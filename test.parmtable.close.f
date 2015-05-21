      PROGRAM testlats

      INCLUDE '/usr/local/lats/LATS-1.1/include/lats.inc'

      DO i=1,10000
        print*, 'Before lats, ',i
        j=latsparmtab('/pcmdi/doutriau/libs/table.lats')
c        print*, i
      ENDDO

      END
