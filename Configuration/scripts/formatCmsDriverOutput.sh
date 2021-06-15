#!/usr/bin/sed -f

# Format the cmsDriver generated config file by applying the following
# - Squeezing Black Lines
# - Replace tabs with spaces
# - Strip trailing whitespaces
# - Insert a space after '#' in lines starting with '#'

# 1 - Squeezing blank lines

# on empty lines, join with next
# Note there is a star in the regexp
:x
/^\n*$/ {
N
bx
}

# now, squeeze all '\n'
s/^\(\n\)*/\1/

# 2 - Replace tabs with spaces

s/\t/    /g

# 3 - Strip trailing whitespaces

s/[[:blank:]]*$//

# 4 - Insert a space after '#' in lines starting with '#'
# Note the M modifier is for multi-line mode

s/^#\(\S\+\)/# \1/m
