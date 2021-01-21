begin
;   load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;;;; setting

grids	= new((/6,2/),"integer")
grids(0,:)	= (/354, 85/)
grids(1,:)	= (/355, 85/)
grids(2,:)	= (/356, 85/)
grids(3,:)	= (/354, 86/)
grids(4,:)	= (/355, 86/)
grids(5,:)	= (/356, 86/)
vars	= (/"PP", "TT", "RH", "WD", "RN", "VP"/)
;;var_names	= (/"Prep", "Temp", "RH", "Wind", "Rn", "Vp"/)
;;lnth	= 17544
lnth	= 8760*3

;;;; file open, write first line

fname	= "/home/ualjxp/NLDAS/ext_yan/ncl_to_forc/sample.txt"
system("rm -f " + fname)
header1	= "6   6   6   6   6   0   6   14   1   0" 
;;header2	= "Time    PP"
hlist	= [/header1/]
write_table(fname, "w", hlist, "%s")
;;write_table(fname, "a", [/header2/], "%s")

;;;;

do i_var=0,2
do i_grid=0,5
   z1	= asciiread("~/NLDAS/ext_yan/x" + grids(i_grid,0) + "y" + grids(i_grid,1) + "z" + vars(i_var) + ".txt", (/lnth,4/), "float")
   z10	= new((/2*lnth,1/), "float")
   z11	= new((/2*lnth,1/), "float")
   do i=0,lnth-1
      z10(2*i,0)	= z1(i,0)
      z10(2*i+1,0)	= z1(i,2)
      z11(2*i,0)	= z1(i,1)
      z11(2*i+1,0)	= z1(i,3)
   end do

   i_grid2	= i_grid + 1
   header2	= vars(i_var) + "  " + i_grid2 + "  " + lnth*2 
   write_table(fname, "a", [/header2/], "%s")
   zlist	= [/z10, z11/]
   write_table(fname, "a", zlist, "%f %f ")
   header3	= "   "
   write_table(fname, "a", [/header3/], "%s")

end do
end do

i_var=3
do i_grid=0,5
   z1   = asciiread("~/NLDAS/ext_yan/x" + grids(i_grid,0) + "y" + grids(i_grid,1) + "z" + vars(i_var) + ".txt", (/lnth,4/), "float")
   z10  = new((/2*lnth,1/), "float")
   z11  = new((/2*lnth,1/), "float")
   do i=0,lnth-1
      z10(2*i,0)        = z1(i,0)
      z10(2*i+1,0)      = z1(i,2)
      z11(2*i,0)        = z1(i,1)
      z11(2*i+1,0)      = z1(i,3)
   end do

   i_grid2      = i_grid + 1
   header2      = vars(i_var) + "  " + i_grid2 + "  " + lnth*2 + " 10"
   write_table(fname, "a", [/header2/], "%s")
   zlist        = [/z10, z11/]
   write_table(fname, "a", zlist, "%f %f ")
   header3      = "   "
   write_table(fname, "a", [/header3/], "%s")
end do

do i_var=4,5
do i_grid=0,5
   z1   = asciiread("~/NLDAS/ext_yan/x" + grids(i_grid,0) + "y" + grids(i_grid,1) + "z" + vars(i_var) + ".txt", (/lnth,4/), "float")
   z10  = new((/2*lnth,1/), "float")
   z11  = new((/2*lnth,1/), "float")
   do i=0,lnth-1
      z10(2*i,0)        = z1(i,0)
      z10(2*i+1,0)      = z1(i,2)
      z11(2*i,0)        = z1(i,1)
      z11(2*i+1,0)      = z1(i,3)
   end do

   i_grid2      = i_grid + 1
   header2      = vars(i_var) + "  " + i_grid2 + "  " + lnth*2
   write_table(fname, "a", [/header2/], "%s")
   zlist        = [/z10, z11/]
   write_table(fname, "a", zlist, "%f %f ")
   header3      = "   "
   write_table(fname, "a", [/header3/], "%s")

end do
end do

;;;;
tstep_veg	= 2*12*3
jday	= (/0, 30, 31, 58, 59, 89, 90, 119, 120, 150, 151, 180, 181, 211, 212, 242, 243, 272, 273, 303, 304, 333, 334, 364, 365, 395, 396, 423, 424, 454, 455, 484, 485, 515, 516, 545, 546, 576, 577, 607, 608, 637, 638, 668, 669, 698, 699, 730, 731, 761, 762, 789, 790, 820, 821, 850, 851, 881, 882, 911, 912, 942, 943, 973, 974, 1003, 1004, 1034, 1035, 1064, 1065, 1095/)
;;jday	= (/0, 30, 31, 58, 59, 89, 90, 119, 120, 150, 151, 180, 181, 211, 212, 242, 243, 272, 273, 303, 304, 333, 334, 365/)
rl_db_nldas	= asciiread("~/NLDAS/ext_yan/ncl_to_forc/rl_nldas.txt",(/14,12/),"float")
lai_db_nldas	= asciiread("~/NLDAS/ext_yan/ncl_to_forc/lai_nldas.txt",(/14,12/),"float")

;; Classes in db_nldas
;; 0: water
;; 1: evergreen needleleaf forest
;; 2: evergreen broadleaf forest
;; 3: deciduous needleleaf forest
;; 4: deciduous broadleaf forest
;; 5: mixed cover
;; 6: woodland
;; 7: wooded grassland
;; 8: closed shrubland
;; 9: open shrubland
;; 10: grassland
;; 11: cropland
;; 12: bare ground
;; 13: urban and built-up

;; Class in NLCD, used in PIHM
;; 0(11): open water <- 0 of nldas
;; 1(21): developed, open space <- 13 of nldas
;; 2(22): developed, low intensity <- 13 of nldas
;; 3(23): developed, medium intensity <- 13 of nldas
;; 4(24): developed, high intensity <- 13 of nldas
;; 5(31): barren land (rock/sand/clay) <- 12 of nldas
;; 6(41): deciduous forest <- 3 of nldas
;; 7(42): evergreen forest <- 1 of nldas
;; 8(43): mixed forest <- 5 of nldas
;; 9(52): shrub/scrub <- 8 of nldas
;; 10(71): grassland/herbaceous <- 10 of nldas
;; 11(81): pasture/hay <- 7 of nldas
;; 12(82): cultivated crops <- 11 of nldas
;; 13(90): woody wetlands <- 6 of nldas
;; 14(95): emergent herbaceous wetlands <- 6 of nldas 

;number_nlcd	= (/11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95/)
number_nlcd	= (/21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95/)
;convert_nldas_nlcd	= (/0, 13, 13, 13, 13, 12, 3, 1, 5, 8, 10, 7, 11, 6, 6/)
convert_nldas_nlcd	= (/13, 13, 13, 13, 12, 3, 1, 5, 8, 10, 7, 11, 6, 6/)
;rl_nlcd	= new((/tstep_veg,15/), "float")
rl_nlcd	= new((/tstep_veg,14/), "float")
;lai_nlcd	= new((/tstep_veg,15/), "float")
lai_nlcd	= new((/tstep_veg,14/), "float")
do ii=0,tstep_veg-1
;do ik=0,14
do ik=0,13
   iii	= round(mod(ii,24)/2 - 0.1,3)
   ;;print(iii)
   ;;iki	= convert_nldas_nlcd(ik)
   ;;print(iki)
   ;;rl_nlcd(ii,ik)	= rl_db_nldas(iki,iii)
   rl_nlcd(ii,ik)	= rl_db_nldas(ik,iii)
   ;;lai_nlcd(ii,ik)	= lai_db_nldas(iki,iii)
   lai_nlcd(ii,ik)	= lai_db_nldas(ik,iii)
end do 
end do

;do i_class_veg=0,14
do i_class_veg=0,13
   i_class_veg2	= i_class_veg + 1
   ;header_veg1	= "LAI  " + number_nlcd(i_class_veg) + "  " + tstep_veg 
   header_veg1	= "LAI  " + i_class_veg2 + "  " + tstep_veg + "  0.0002" 
   write_table(fname, "a", [/header_veg1/], "%s")
   lai_temp1	= jday
   lai_temp2	= lai_nlcd(:,i_class_veg)
   zlist_veg	= [/lai_temp1, lai_temp2/]
   write_table(fname, "a", zlist_veg, "%d%8.4f")
   header_veg2	= "   "
   write_table(fname, "a", [/header_veg2/], "%s")
end do

;do i_class_veg=0,14
do i_class_veg=0,13
   i_class_veg2	= i_class_veg + 1
   ;header_veg1	= "RL  " + number_nlcd(i_class_veg) + "  " + tstep_veg
   header_veg1	= "RL  " + i_class_veg2 + "  " + tstep_veg
   write_table(fname, "a", [/header_veg1/], "%s")
   rl_temp1	= jday
   rl_temp2	= rl_nlcd(:,i_class_veg)
   zlist_veg	= [/rl_temp1, rl_temp2/]
   write_table(fname, "a", zlist_veg, "%d%8.4f")
   header_veg2	= "   "
   write_table(fname, "a", [/header_veg2/], "%s")
end do

;;;;

header_MF	= "MF  1  2"
MF1	= "0   0.0"
MF2	= "365   0.0"
write_table(fname, "a", [/header_MF/], "%s")
write_table(fname, "a", [/MF1/], "%s")
write_table(fname, "a", [/MF2/], "%s")

;;;;

end
