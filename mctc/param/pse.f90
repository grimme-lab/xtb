module mctcpar_pse
   implicit none
   character(len=2), parameter :: pse(118) = [ &
      & 'h ','he', &
      & 'li','be','b ','c ','n ','o ','f ','ne', &
      & 'na','mg','al','si','p ','s ','cl','ar', &
      & 'k ','ca', &
      & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
      &           'ga','ge','as','se','br','kr', &
      & 'rb','sr', &
      & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
      &           'in','sn','sb','te','i ','xe', &
      & 'cs','ba','la', &
      & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
      & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
      &           'tl','pb','bi','po','at','rn', &
      & 'fr','ra','ac', &
      & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
      & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
      &           'nh','fl','mc','lv','ts','og' ]

   integer, parameter :: metal(94) = [&
      &  0,                                                 0, &! H-He
      &  1, 1,                               1, 0, 0, 0, 0, 0, &! Li-Ne
      &  1, 1,                               1, 0, 0, 0, 0, 0, &! Na-Ar
      &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! K-Kr
      &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Rb-Xe
      &  1, 1, &! Cs/Ba
      &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &!La-Lu
      &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Lu-Rn
      &  1, 1, 1, 1, 1, 1, 1, 1  ]! Fr-Pu

end module mctcpar_pse
