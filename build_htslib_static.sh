#!/bin/bash
install_prefix=""

## must be in htslib source directory
# Ensure Makefile and sam.c exist in the current directory
if [[ ! -f Makefile || ! -f sam.c ]]; then
  echo "Error: this script must be run in the htslib source directory." >&2
  exit 1
fi

# Parse optional install_prefix argument
if [ "$#" -ge 1 ]; then
  install_prefix=$1
fi

# Check if install-lib-static target exists in the Makefile
if ! grep -q "^install-lib-static:" Makefile; then
  echo "Adding install-lib-static target to the Makefile..."
  perl -i -pe 's/^install:/install-lib-static: lib-static \$(BUILT_PROGRAMS)\
\t\$(INSTALL_DIR) \$(DESTDIR)\$(bindir) \$(DESTDIR)\$(includedir) \$(DESTDIR)\$(includedir)\/htslib \$(DESTDIR)\$(libdir)\
\t\$(INSTALL_PROGRAM) \$(BUILT_PROGRAMS) \$(DESTDIR)\$(bindir)\
\t\$(INSTALL_DATA) libhts.a \$(DESTDIR)\$(libdir)\/libhts.a\
\t\$(INSTALL_DATA) \$(SRC)htslib\/\*.h \$(DESTDIR)\$(includedir)\/htslib\
\ninstall:/' Makefile
fi

# Build the static library
echo "Building the static library..."
make -j6 lib-static

# Install the static library and headers if install_prefix is provided
if [ -n "$install_prefix" ]; then
  echo "Installing static library and headers to prefix: $install_prefix..."
  make install-lib-static prefix="$install_prefix"
else
  echo "Installation prefix not provided. Build complete."
fi


#perl -0777 -pi -e 's/(^install:)/install-lib-static: lib-static installdirs\n\t$(INSTALL_DATA)\
# libhts.a $(DESTDIR)$(libdir)\/libhts.a\n\t$(INSTALL_DATA) $(SRC)htslib\/\*.h\
#  $(DESTDIR)$(includedir)\/htslib\n\n\1/' Makefile
