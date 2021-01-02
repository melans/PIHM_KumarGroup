#echo "$PWD" > /glade/work/elansary/somerville_150/5055/12h/rslts.csv
##echo "`ls -al`" >> rslts.csv
##echo 'Time,Points,R2,Sprmn2,NSE,mNSE,NSE_cfs,mNSE_cfs,Setup' >> /glade/work/elansary/somerville_150/5055/12h/rslts.csv;
#echo -e 'Ttime\tPoints\tR2\tSprmn2\tNSE\tmNSE\tNSE_cfs\tmNSE_cfs\tSetup' >> /glade/work/elansary/somerville_150/5055/12h/rslts.csv;
#files=(`wc -l */*x1|awk '$1>365*2&&$2!="total"{split($2,f,"/");printf f[1]" "}'`);

echo "$PWD" > rslts.csv
echo -e 'Ttime\tPoints\tR2\tSprmn2\tNSE\tmNSE\tNSE_cfs\tmNSE_cfs\tSetup' >> rslts.csv;
files=*x1

for f in "${files[@]}";do 
	Rscript stats.r `realpath $f` --vanilla #>/dev/null 2>/dev/null;
done;
cat rslts.csv;

#mail -s "PIHM: $PWD (`date`)" -a rslts.csv -a stats.r -a stats.sh melans@tamu.edu < .
