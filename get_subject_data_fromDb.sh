#!/bin/bash
read -r -d '' copysel << EOT
 copy (select 'P'||trim(to_char(s.id, '00000')) as subj_id, brnum as "BrNum", dx as "Dx", sex as "Sex", 
        race as "Race", age as "Age", coalesce(pmi, 0.0) as "PMI", coalesce(mod, '.') as "MoD", dropped 
        from subjects s, dx d where d.id=s.dx_id ORDER BY 1) 
    to stdout with (format 'csv',header true)
EOT
psql -h srv16 -d rse -c "$copysel" | csv2tab 
