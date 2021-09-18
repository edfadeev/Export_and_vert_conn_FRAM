### primer clipping from 16S libraries produced using V4V5 primers
###Author: Eddie Fadeev

# set variables, program and script locations
FILELOCATION="/scratch2/efadeev/PS99/dada2"
NSAMPLE="72"

# scripts
CLIP_qsub="/scratch2/efadeev/aphros/ampliconNGS/clipping_aphros.sh"

###step 0: shorten file names
# save original file names
ls -1v ./Original/*R1_001.fastq > ./originalR1
ls -1v ./Original/*R2_001.fastq > ./originalR2

mkdir Renamed

# copy original files to new file names
new=1
for i in $(ls -1v ./Original/*R1_001.fastq)
do
cp ${i} ./Renamed/${new}"_R1.fastq"
((new++))
done

new=1
for i in $(ls -1v ./Original/*R2_001.fastq)
do
cp ${i} ./Renamed/${new}"_R2.fastq"
((new++))
done

# check that the renaming was done correctly
ls -1v ./Renamed/[0-9]*_R1.fastq > ./renamedR1
ls -1v ./Renamed/[0-9]*_R2.fastq > ./renamedR2

paste ./originalR1 ./renamedR1 > ./fileID_R1
paste ./originalR2 ./renamedR2 > ./fileID_R2

#the following commands schould not give any output
while read line ; do
diff $(echo "$line")
done < ./fileID_R1

while read line ; do
diff $(echo "$line")
done < ./fileID_R2


# create directory for qsub output files
mkdir Logfiles

###step 1: primer clipping 

# bacterial primer V4-V5
FBC=^GTGYCAGCMGCCGCGGTAA # forward primer
RBC=^CCGYCAATTYMTTTRAGTTT # reverse primer
OFWD=18 # length of forward primer (19) - 1
OREV=19 # length of reverse primer (20) - 1
ERROR=0.16


mkdir Clipped

qsub -d ${FILELOCATION} -t 1-${NSAMPLE} -o ${FILELOCATION}/Logfiles -e ${FILELOCATION}/Logfiles -v FBC=${FBC},RBC=${RBC},OFWD=${OFWD},OREV=${OREV},ERROR=${ERROR} ${CLIP_qsub}

# cleaning up directories
mkdir ./Clipped/Clipped_logs
mv ./Clipped/*.log ./Clipped/Clipped_logs/
  mv ./Clipped/*.info ./Clipped/Clipped_logs/
  mkdir ./Clipped/rf
mv Clipped/*rf* Clipped/rf/
  
  #generate sample names
  cd Clipped
ls -1v *.fastq > ../sample_names.txt
cd ..