#! /usr/bin/awk -f

/Datos extraidos del fichero:/ { file = $NF }

/Bond critical points/,/^ *$/ {
   if ( ($2=="P" && $4=="N") || ($2=="N" && $4=="P") ) {
      dist = $6 + $7
      rho = $11
      del2 = $12
      kinG = $13
      lam1 = $14
      lam2 = $15
      lam3 = $16
      printf "%-45s", file
      print dist, rho, del2, kinG, lam1, lam2, lam3
      }
   }
