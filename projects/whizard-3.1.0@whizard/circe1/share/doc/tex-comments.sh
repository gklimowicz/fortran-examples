#! /usr/bin/awk -f
# tex-comments.sh --

/^@begin docs / { code = 0 }
/^@begin code / { code = 1 }

code && /^@text .*![:$]/ {
  if (match($0, /!:.*$/)) {
    printf("%s\n", substr($0, 1, RSTART-1)) 
    printf("@literal ! {\\setupmodname %s}\n", substr($0, RSTART+2))
    next
  }
  if (match($0, /!\$.*$/)) {
    printf("%s\n", substr($0, 1, RSTART-1)) 
    printf("@literal ! {\\setupmodname$ %s $}\n", substr($0, RSTART+2))
    next
  }
}

# Hide a trick for Poor Man's Elemental Procedures
code { gsub(/`'_/, "_") }
  
{ print }
