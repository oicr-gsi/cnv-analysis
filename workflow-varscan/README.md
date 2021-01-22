# varscan

Varscan 2.2, workflow for calling SNVs and CVs

## Overview

Varscan workflow, calls Copy Number Variants and SNVs
Creation of mpileups and calling variants are done with parallel processing

![varscan outputs](docs/Screenshot_Varscan.png)

## Dependencies

* [picard 2.21.2](https://broadinstitute.github.io/picard)
* [varscan 2.4.2](http://varscan.sourceforge.net)
* [samtools 0.1.19](http://www.htslib.org/)
* [rstats 3.6](http://cran.utstat.utoronto.ca/src/base/R-3/R-3.6.1.tar.gz)


## Usage

### Cromwell
```
java -jar cromwell.jar run varscan.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputTumor`|File|input .bam file for tumor sample
`inputNormal`|File|input .bam file for normal sample
`inputTumorIndex`|File|input .bai file for tumor sample
`inputNormalIndex`|File|input .bai file for normal sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|""|Output file(s) prefix
`bedIntervalsPath`|String|""|Path to a .bed file used for targeted variant calling
`chromRegions`|Array[String]|["chr1:1-249250621", "chr2:1-243199373", "chr3:1-198022430", "chr4:1-191154276", "chr5:1-180915260", "chr6:1-171115067", "chr7:1-159138663", "chr8:1-146364022", "chr9:1-141213431", "chr10:1-135534747", "chr11:1-135006516", "chr12:1-133851895", "chr13:1-115169878", "chr14:1-107349540", "chr15:1-102531392", "chr16:1-90354753", "chr17:1-81195210", "chr18:1-78077248", "chr19:1-59128983", "chr20:1-63025520", "chr21:1-48129895", "chr22:1-51304566", "chrX:1-155270560", "chrY:1-59373566", "chrM:1-16571"]|Regions used for scattering tasks, need to be assembly-specific


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`expandRegions.jobMemory`|Int|4|Memory for this task in GB
`makePileups.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Reference fasta file, path depends on the respective module
`makePileups.modules`|String|"samtools/0.1.19 hg19/p13"|required modules
`makePileups.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|path to samtools
`makePileups.jobMemory`|Int|18|memory for this job, in Gb
`makePileups.timeout`|Int|40|Timeout in hours, needed to override imposed limits
`runVarscanCNV.pValue`|Float|0.05|p-value for cnv calling, default is 0.05
`runVarscanCNV.jobMemory`|Int|20|Memory in Gb for this job
`runVarscanCNV.javaMemory`|Int|6|Memory in Gb for Java
`runVarscanCNV.logFile`|String|"VARSCAN_CNV.log"|File for logging Varscan messages
`runVarscanCNV.varScan`|String|"$VARSCAN_ROOT/VarScan.jar"|path to varscan .jar file
`runVarscanCNV.modules`|String|"varscan/2.4.2 java/8"|Names and versions of modules
`runVarscanCNV.timeout`|Int|40|Timeout in hours, needed to override imposed limits
`getSnvNative.pValue`|Float|0.05|somatic p-value for SNV calling, default is 0.05
`getSnvNative.jobMemory`|Int|20|Memory in Gb for this job
`getSnvNative.javaMemory`|Int|6|Memory in Gb for Java
`getSnvNative.minCoverage`|Int|8|Minimum coverage in normal and tumor to call variant [8]
`getSnvNative.minCoverageNormal`|Int|8|Minimum coverage in normal to call somatic [8]
`getSnvNative.minCoverageTumor`|Int|6|Minimum coverage in tumor to call somatic [6]
`getSnvNative.minVarFreq`|Float|0.1|Minimum variant frequency to call a heterozygote [0.10]
`getSnvNative.minFreqForHom`|Float|0.75|Minimum frequency to call homozygote [0.75]
`getSnvNative.normalPurity`|Float|1.0|Estimated purity (non-tumor content) of normal sample [1.00]
`getSnvNative.tumorPurity`|Float|1.0|Estimated purity (tumor content) of normal sample [1.00]
`getSnvNative.pValueHet`|Float|0.99|p-value threshold to call a heterozygote [0.99]
`getSnvNative.strandFilter`|Int|0|If set to 1, removes variants with >90% strand bias
`getSnvNative.validation`|Int|0|If set to 1, outputs all compared positions even if non-variant
`getSnvNative.outputVcf`|Int|0|Flag that when set to 1 indicates that we need results in vcf format
`getSnvNative.logFile`|String|"VARSCAN_SNV.log"|File for logging Varscan messages
`getSnvNative.varScan`|String|"$VARSCAN_ROOT/VarScan.jar"|path to varscan .jar file
`getSnvNative.modules`|String|"varscan/2.4.2 java/8"|Names and versions of modules
`getSnvNative.timeout`|Int|40|Timeout in hours, needed to override imposed limits
`getSnvVcf.pValue`|Float|0.05|somatic p-value for SNV calling, default is 0.05
`getSnvVcf.jobMemory`|Int|20|Memory in Gb for this job
`getSnvVcf.javaMemory`|Int|6|Memory in Gb for Java
`getSnvVcf.minCoverage`|Int|8|Minimum coverage in normal and tumor to call variant [8]
`getSnvVcf.minCoverageNormal`|Int|8|Minimum coverage in normal to call somatic [8]
`getSnvVcf.minCoverageTumor`|Int|6|Minimum coverage in tumor to call somatic [6]
`getSnvVcf.minVarFreq`|Float|0.1|Minimum variant frequency to call a heterozygote [0.10]
`getSnvVcf.minFreqForHom`|Float|0.75|Minimum frequency to call homozygote [0.75]
`getSnvVcf.normalPurity`|Float|1.0|Estimated purity (non-tumor content) of normal sample [1.00]
`getSnvVcf.tumorPurity`|Float|1.0|Estimated purity (tumor content) of normal sample [1.00]
`getSnvVcf.pValueHet`|Float|0.99|p-value threshold to call a heterozygote [0.99]
`getSnvVcf.strandFilter`|Int|0|If set to 1, removes variants with >90% strand bias
`getSnvVcf.validation`|Int|0|If set to 1, outputs all compared positions even if non-variant
`getSnvVcf.logFile`|String|"VARSCAN_SNV.log"|File for logging Varscan messages
`getSnvVcf.varScan`|String|"$VARSCAN_ROOT/VarScan.jar"|path to varscan .jar file
`getSnvVcf.modules`|String|"varscan/2.4.2 java/8"|Names and versions of modules
`getSnvVcf.timeout`|Int|40|Timeout in hours, needed to override imposed limits
`mergeCNV.jobMemory`|Int|6|memory in GB for this job
`mergeCNV.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`mergeSNP.jobMemory`|Int|6|memory in GB for this job
`mergeSNP.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`mergeIND.jobMemory`|Int|6|memory in GB for this job
`mergeIND.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`mergeSNPvcf.modules`|String|"picard/2.21.2 hg19/p13"|modules needed for this task
`mergeSNPvcf.seqDictionary`|String|"$HG19_ROOT/hg19_random.dict"|.dict file for the reference in use
`mergeSNPvcf.jobMemory`|Int|12|memory in GB for this job
`mergeSNPvcf.javaMemory`|Int|8|memory in GB for java VM
`mergeSNPvcf.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`mergeINDvcf.modules`|String|"picard/2.21.2 hg19/p13"|modules needed for this task
`mergeINDvcf.seqDictionary`|String|"$HG19_ROOT/hg19_random.dict"|.dict file for the reference in use
`mergeINDvcf.jobMemory`|Int|12|memory in GB for this job
`mergeINDvcf.javaMemory`|Int|8|memory in GB for java VM
`mergeINDvcf.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`smoothData.varScan`|String|"$VARSCAN_ROOT/VarScan.jar"|Path to VarScan jar file
`smoothData.modules`|String|"varscan/2.4.2 java/8 rstats/3.6"|Modules for this job
`smoothData.min_coverage`|Int|20|Fine-tuning parameter for VarScan
`smoothData.max_homdel_coverage`|Int|5|Max coverage form homozygous deletion, default is 5
`smoothData.min_tumor_coverage`|Int|10|Min coverage in tumor sample, default is 10
`smoothData.del_threshold`|Float|0.25|Fine-tuning parameter for VarScan
`smoothData.amp_threshold`|Float|0.25|Amplification threshold to report, default is 0.25
`smoothData.min_region_size`|Int|10|Fine-tuning parameter for VarScan
`smoothData.recenter_up`|Int|0|Fine-tuning parameter for VarScan
`smoothData.recenter_down`|Int|0|Fine-tuning parameter for VarScan
`smoothData.jobMemory`|Int|16|Memory in Gb for this job
`smoothData.javaMemory`|Int|6|Memory in Gb for Java


### Outputs

Output | Type | Description
---|---|---
`resultCnvFile`|File?|file with CNV calls, smoothed
`resultSnpFile`|File?|file with SNPs, native varscan format
`resultIndelFile`|File?|file with Indel calls, native varscan format
`resultSnpVcfFile`|File?|file with SNPs, vcf format
`resultIndelVcfFile`|File?|file with Indels, vcf format


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
