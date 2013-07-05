#! /usr/bin/gawk -f

BEGIN{ g = 0 }

/bravais-lattice index/ { ibrav = $4 }
/lattice parameter \(alat\)/ { alat = $5 }
/crystal axes:/ {
   g++;
   getline; a[g,1,1]=$4; a[g,1,2]=$5; a[g,1,3]=$6;
   getline; a[g,2,1]=$4; a[g,2,2]=$5; a[g,2,3]=$6;
   getline; a[g,3,1]=$4; a[g,3,2]=$5; a[g,3,3]=$6;
   }

END {
   file = "tmppw2cell"
   printf "%d\n", ibrav >  file
   printf "%s\n", alat  >> file
   for (i=1; i<=3; i++) {
      printf "%s %s %s\n", a[g,i,1], a[g,i,2], a[g,i,3]  >> file
      }
   order = sprintf("cr_pw2cell.m %s", file)
   system(order)
   }
