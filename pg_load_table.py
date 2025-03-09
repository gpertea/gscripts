#!/bin/env python
import pandas as pd
import psycopg2
import os
import sys
import getopt
import csv

def csv_to_pg_table(csv_file, table_name, conn_params):
    """Loads a CSV/TSV file into PostgreSQL, inferring types and creating the table."""
    conn = None
    cur = None
    try:
        # Determine delimiter by reading the first line
        with open(csv_file, 'r') as f:
            first_line = f.readline()
            delimiter = '\t' if '\t' in first_line else ','
            has_header = True  # Assume the file has a header

        # Use pandas to read the file
        df = pd.read_csv(csv_file, sep=delimiter, engine='python', header=0 if has_header else None)

        conn = psycopg2.connect(**conn_params)
        cur = conn.cursor()

        # Check if table exists
        cur.execute(f"SELECT EXISTS (SELECT 1 FROM pg_tables WHERE schemaname = 'public' AND tablename = '{table_name}');")
        table_exists = cur.fetchone()[0]

        if not table_exists:
            # Create table with inferred column types
            create_table_sql = f"CREATE TABLE {table_name} ("
            column_names = df.columns if has_header else [f"col{i+1}" for i in range(df.shape[1])]
            for col_name, col_type in zip(column_names, df.dtypes):
                if pd.api.types.is_integer_dtype(col_type):
                    sql_type = "INTEGER"
                elif pd.api.types.is_float_dtype(col_type):
                    sql_type = "REAL"
                elif pd.api.types.is_datetime64_any_dtype(col_type):
                    sql_type = "TIMESTAMP"
                elif str(col_type) == 'bool':
                    sql_type = "BOOLEAN"
                else:
                    sql_type = "TEXT"
                create_table_sql += f"\"{col_name}\" {sql_type}, "
            create_table_sql = create_table_sql[:-2] + ")"
            cur.execute(create_table_sql)
            conn.commit()

        # Load data using COPY
        with open(csv_file, 'r') as f:
            cur.copy_expert(f"COPY {table_name} FROM STDIN WITH (FORMAT CSV, HEADER {has_header}, DELIMITER E'{delimiter}')", f)
            conn.commit()
            print(f"Data from '{csv_file}' loaded into '{table_name}'.")

    except (Exception, psycopg2.Error) as error:
        print(f"Error: {error}")
        if conn:
            conn.rollback()
    finally:
        if cur:
            cur.close()
        if conn:
            conn.close()

def get_conn_params(args):
    """Gets connection parameters from environment, .pgpass, and command line."""

    conn_params = {}

    conn_params['host'] = os.environ.get('PGHOST')
    conn_params['dbname'] = os.environ.get('PGDATABASE')
    conn_params['user'] = os.environ.get('PGUSER')
    conn_params['port'] = os.environ.get('PGPORT')
    conn_params['password'] = os.environ.get('PGPASSWORD')

    pgpass_file = os.path.expanduser("~/.pgpass")
    if os.path.exists(pgpass_file):
        with open(pgpass_file, "r") as f:
            for line in f:
                parts = line.strip().split(":")
                if len(parts) == 5:
                    pg_host = conn_params.get('host', 'localhost')
                    pg_user = conn_params.get('user')
                    pg_dbname = conn_params.get('dbname')
                    pg_port = conn_params.get('port', '5432')
                    if (parts[0] == pg_host or parts[0] == '*') and \
                       (parts[1] == pg_port or parts[1] == '*') and \
                       (parts[2] == pg_dbname or parts[2] == '*') and \
                       (parts[3] == pg_user or parts[3] == '*'):
                           conn_params['host'] = parts[0] if parts[0] != '*' else pg_host
                           conn_params['port'] = parts[1] if parts[1] != '*' else pg_port
                           conn_params['dbname'] = parts[2] if parts[2] != '*' else pg_dbname
                           conn_params['user'] = parts[3] if parts[3] != '*' else pg_user
                           conn_params['password'] = parts[4]
                           break

    conn_params.update(args)

    if 'port' in conn_params and conn_params['port'] is not None:
        try:
            conn_params['port'] = int(conn_params['port'])
        except ValueError:
            print("Error: Port must be an integer.")
            sys.exit(1)

    return conn_params


if __name__ == "__main__":
    csv_file = None
    table_name = None
    conn_params = {}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:d:U:p:W:",
                                   ["host=", "database=", "user=", "port=", "password="])
        if len(args) == 2:
            csv_file = args[0]
            table_name = args[1]
        else:
            print("Usage: pg_load_table.py csv_file table_name [-h host] [-d database] [-U user] [-p port] [-W password]")
            sys.exit(1)

        args_dict = {}
        for opt, arg in opts:
            if opt in ("-h", "--host"):
                args_dict['host'] = arg
            elif opt in ("-d", "--database"):
                args_dict['dbname'] = arg
            elif opt in ("-U", "--user"):
                args_dict['user'] = arg
            elif opt in ("-p", "--port"):
                args_dict['port'] = arg
            elif opt in ("-W", "--password"):
                args_dict['password'] = arg

    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    conn_params = get_conn_params(args_dict)
    csv_to_pg_table(csv_file, table_name, conn_params)
