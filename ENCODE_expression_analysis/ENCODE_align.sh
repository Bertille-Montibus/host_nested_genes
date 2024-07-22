# Author: James Cain 
# Downloading ENCODE data 


# Download files from ENCODE Cart
xargs -L 1 curl -O -J -L < files.txt


# Provide first and second read here in array format 
array=(ENCFF576OBS
ENCFF541KUW
ENCFF338DKW
ENCFF540SNP
ENCFF850ZLY
ENCFF456MMS
ENCFF016TGP
ENCFF034BHU
ENCFF402BWO
ENCFF187OKV
ENCFF650JAM
ENCFF002BEB
ENCFF002BEF
ENCFF482WYS
ENCFF926YPC
ENCFF734ZAD
ENCFF322RPT
ENCFF433PKC
ENCFF911WYB
ENCFF912LIX
ENCFF044VER
ENCFF028DUO
ENCFF709FHN
ENCFF419GVS
ENCFF170RHF
ENCFF592VVB
ENCFF770NYA
ENCFF464TEM)

array2=( ENCFF063OPQ
ENCFF386MOY
ENCFF721RGF
ENCFF999PRA
ENCFF897IUQ
ENCFF716WNR
ENCFF604DIX
ENCFF788AET
ENCFF450JES
ENCFF283RUU
ENCFF803DXA
ENCFF002BEC
ENCFF002BEG
ENCFF058MGQ
ENCFF111IRS
ENCFF261RWK
ENCFF782AHJ
ENCFF546XEW
ENCFF723QXK
ENCFF057RDO
ENCFF640PYL
ENCFF470RWW
ENCFF681HNP
ENCFF135CVY
ENCFF437XFH
ENCFF359HIQ
ENCFF076IRZ
ENCFF221QNJ )

# Generation of FastQC Reports
mkdir raw
mkdir FastQC
for ID in $array; do
	mv ${ID}.fastq.gz raw
	fastqc raw/${ID}.fastq.gz --outdir FastQC/
	done
for ID in $array2; do
  mv ${ID}.fastq.gz raw
  fastqc raw/${ID}.fastq.gz --outdir FastQC/
  done

# Kallisto index build and alignment + counting 
kallisto index -i gencode.v38.alltranscripts.idx gencode.v38.transcripts.fa.gz
mkdir kalliso_output_gencode
for ((i=0;i<${#array[@]}; i++)); do
    mkdir kalliso_output_gencode/${array[$i]}
    kallisto quant -i homo_sapiens_kallisto_index/gencode.v38.alltranscripts.idx -o kalliso_output_gencode/${array[$i]} raw/${array[$i]}.fastq.gz raw/${array2[$i]}.fastq.gz
done

# tspex use to generate tau 
ts=(tau
 gini
 simpson
 shannon_specificity
 roku_specificity
 spm_dpm
 js_specificity_dpm)
for t in $ts; do
tspex normalisedcounts_gencode.csv tissue_specific_gencode/${t}.tsv ${t}
done


## ENCODE ANALYSIS WITH SIMILAR MOUSE DATA ## 
# Repeating Alignment with Mouse Data

mkdir mouse_raw
xargs -L 1 curl -O -J -L < Mouse_ENCODE.txt

# Generation of mouse arrays for paired end reads
# Provide first and second read here:
array=( ENCFF563FDS
ENCFF307YNT
ENCFF001IYI
ENCFF001IYK
ENCFF001IYJ
ENCFF001IYL
ENCFF001JTI
ENCFF001JTM
ENCFF161LEK
ENCFF492PRP
ENCFF848QCK
ENCFF432KKN
ENCFF786ZKB
ENCFF682XSC
ENCFF871XHK
ENCFF104UFH
ENCFF547BPL
ENCFF128KGA
ENCFF463WEH
ENCFF286WTQ
ENCFF445AWP end)
array2=( ENCFF470ULL
ENCFF731TMT
ENCFF001IYM
ENCFF001IYO
ENCFF001IYN
ENCFF001IYP
ENCFF001JTW
ENCFF001JWN
ENCFF516HOO
ENCFF581OEV
ENCFF511PCY
ENCFF859JTH
ENCFF517RDO
ENCFF690HKC
ENCFF952JKH
ENCFF126WYO
ENCFF554BIM
ENCFF510DLJ
ENCFF312OKA
ENCFF358NPU
ENCFF958CHE end)

kallisto index -i gencode.vM27.transcripts.idx gencode.vM27.transcripts.fa.gz
mkdir mouse_kalliso_output_gencode
for ((i=0;i<${#array[@]}; i++)); do
    mkdir mouse_kalliso_output_gencode/${array[$i]}
    kallisto quant -i mus_musculus_kallisto_index/gencode.vM27.transcripts.idx -o mouse_kalliso_output_gencode/${array[$i]} mouse_raw/${array[$i]}.fastq.gz mouse_raw/${array2[$i]}.fastq.gz
done

# tspex use for mouse data 
ts=(tau
 gini
 simpson
 shannon_specificity
 roku_specificity
 spm_dpm
 js_specificity_dpm)
mkdir tissue_specific_mouse_gencode
for t in $ts; do
tspex mouse_normalisedcounts_gencode.csv tissue_specific_mouse_gencode/${t}.tsv ${t}
done
