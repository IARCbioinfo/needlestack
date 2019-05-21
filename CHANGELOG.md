# Change Log

## [v1.1](https://github.com/IARCbioinfo/needlestack/tree/v1.1) (2019-05-21)
[Full Changelog](https://github.com/IARCbioinfo/needlestack/compare/v1.0...v1.1)

**Implemented enhancements:**

- add needlestack logo in log.info [\#182](https://github.com/IARCbioinfo/needlestack/issues/182)
- add INFO field containing allelic frequency [\#174](https://github.com/IARCbioinfo/needlestack/issues/174)
- Add an option to only plot somatic mutations [\#149](https://github.com/IARCbioinfo/needlestack/issues/149)
- Make the main 3 processes pipe-friendly to avoid intermediate outputs [\#147](https://github.com/IARCbioinfo/needlestack/issues/147)
- Re-write pileup2baseindel.pl in C++ [\#146](https://github.com/IARCbioinfo/needlestack/issues/146)
- Flag possible cross-sample contamination of normal DNA in tumor samples [\#138](https://github.com/IARCbioinfo/needlestack/issues/138)
- Put ./. genotype when there is no power to identify a variant [\#137](https://github.com/IARCbioinfo/needlestack/issues/137)
- Implement Tumor-Normal pair somatic variant calling  [\#133](https://github.com/IARCbioinfo/needlestack/issues/133)
- Add alignment plot in the PDF [\#73](https://github.com/IARCbioinfo/needlestack/issues/73)
- Make a script that would run needlestack without nextflow [\#72](https://github.com/IARCbioinfo/needlestack/issues/72)

**Fixed bugs:**

- max\_dp default value in readme is wrong [\#176](https://github.com/IARCbioinfo/needlestack/issues/176)
- --help should exit 0 [\#161](https://github.com/IARCbioinfo/needlestack/issues/161)
- Manage bed specification for both bed and region [\#157](https://github.com/IARCbioinfo/needlestack/issues/157)
- Using  renamed symlinks as input BAM files doesn't work with --use\_file\_name option [\#152](https://github.com/IARCbioinfo/needlestack/issues/152)

**Closed issues:**

- How to identify each sample? [\#166](https://github.com/IARCbioinfo/needlestack/issues/166)
- Github readme != Docker readme [\#159](https://github.com/IARCbioinfo/needlestack/issues/159)

## [v1.0](https://github.com/IARCbioinfo/needlestack/tree/v1.0) (2016-08-03)
[Full Changelog](https://github.com/IARCbioinfo/needlestack/compare/v0.3...v1.0)

**Implemented enhancements:**

- Manage the three possible genotypes in vcf [\#130](https://github.com/IARCbioinfo/needlestack/issues/130)
- The graph showing AF vs log10\(qval\) should show phred-scaled qvalues [\#121](https://github.com/IARCbioinfo/needlestack/issues/121)

**Fixed bugs:**

- Contours seem to be incorrect [\#128](https://github.com/IARCbioinfo/needlestack/issues/128)
- correct file name extraction for sample name [\#126](https://github.com/IARCbioinfo/needlestack/issues/126)
- Let min\_qval be equal to 0 [\#119](https://github.com/IARCbioinfo/needlestack/issues/119)
- plot improved error rate confidence interval  [\#117](https://github.com/IARCbioinfo/needlestack/issues/117)

**Closed issues:**

- QUAL should not be reported as Inf in VCF when q-value=0 [\#125](https://github.com/IARCbioinfo/needlestack/issues/125)
- Add pipeline execution DAG in README [\#123](https://github.com/IARCbioinfo/needlestack/issues/123)

## [v0.3](https://github.com/IARCbioinfo/needlestack/tree/v0.3) (2016-05-03)
[Full Changelog](https://github.com/IARCbioinfo/needlestack/compare/v0.2...v0.3)

**Implemented enhancements:**

- color points by qvalues in regression plot [\#85](https://github.com/IARCbioinfo/needlestack/issues/85)
- Add an option to directly input a region for calling in the command line [\#71](https://github.com/IARCbioinfo/needlestack/issues/71)
- Improve the bed split method [\#47](https://github.com/IARCbioinfo/needlestack/issues/47)
- Change the number of entry in the INFO and FORMAT VCF fields [\#108](https://github.com/IARCbioinfo/needlestack/issues/108)
- Add contour lines for a set of qvalues in the plot [\#100](https://github.com/IARCbioinfo/needlestack/issues/100)
- Add an option to choose output VCF file name \(--out\_vcf?\) [\#81](https://github.com/IARCbioinfo/needlestack/issues/81)
- Change the way we publish new version [\#69](https://github.com/IARCbioinfo/needlestack/issues/69)
- Make the stable docker file more stable [\#68](https://github.com/IARCbioinfo/needlestack/issues/68)
- Add more tests in CircleCI [\#55](https://github.com/IARCbioinfo/needlestack/issues/55)
- Remove unnecessary intermediate outputs [\#51](https://github.com/IARCbioinfo/needlestack/issues/51)
- In the absence of a bed file the pipeline should run on the full reference genome [\#39](https://github.com/IARCbioinfo/needlestack/issues/39)
- Improve R script command line parsing [\#38](https://github.com/IARCbioinfo/needlestack/issues/38)
- Add version numbers in VCF output [\#20](https://github.com/IARCbioinfo/needlestack/issues/20)

**Fixed bugs:**

- VCF files have to be sorted [\#110](https://github.com/IARCbioinfo/needlestack/issues/110)
- Sometimes large number in VCF files are written in scientific notations [\#109](https://github.com/IARCbioinfo/needlestack/issues/109)
- error when coverage is null for every bam file [\#99](https://github.com/IARCbioinfo/needlestack/issues/99)
- Calling doesn't work when a region contains only T in the reference [\#96](https://github.com/IARCbioinfo/needlestack/issues/96)
- Check that BAM folder contains bam files  [\#66](https://github.com/IARCbioinfo/needlestack/issues/66)
- Check if the gzi is present if the ref is gz [\#65](https://github.com/IARCbioinfo/needlestack/issues/65)
- Verify the user inputs are correct [\#42](https://github.com/IARCbioinfo/needlestack/issues/42)

## [v0.2](https://github.com/IARCbioinfo/needlestack/tree/v0.2) (2015-10-19)
[Full Changelog](https://github.com/IARCbioinfo/needlestack/compare/v0.1...v0.2)

**Implemented enhancements:**

- Add logo image [\#62](https://github.com/IARCbioinfo/needlestack/issues/62)
- add --no\_indel option [\#56](https://github.com/IARCbioinfo/needlestack/issues/56)
- Correct english typos in readme, help and log [\#53](https://github.com/IARCbioinfo/needlestack/issues/53)
- The pipeline randomly crashes with java.nio.file.NoSuchFileException: XXX\_empty.pdf [\#49](https://github.com/IARCbioinfo/needlestack/issues/49)
- Add information about the pipeline in the log [\#41](https://github.com/IARCbioinfo/needlestack/issues/41)
- Add program usage when launched with --help [\#40](https://github.com/IARCbioinfo/needlestack/issues/40)
- Change the way chromosome length is calculated [\#34](https://github.com/IARCbioinfo/needlestack/issues/34)
- Option `all\_sites` should rather be called `all\_SNVs` [\#33](https://github.com/IARCbioinfo/needlestack/issues/33)
- Choose a better name for the pipeline and change file names accordingly [\#31](https://github.com/IARCbioinfo/needlestack/issues/31)
- Add contigs in VCF header [\#25](https://github.com/IARCbioinfo/needlestack/issues/25)
- Move to IARC-bioinfo organisation repo [\#22](https://github.com/IARCbioinfo/needlestack/issues/22)
- Add a zoomed regression plot [\#21](https://github.com/IARCbioinfo/needlestack/issues/21)

**Fixed bugs:**

- The pipeline randomly crashes with java.nio.file.NoSuchFileException: XXX\\_empty.pdf [\#49](https://github.com/IARCbioinfo/needlestack/issues/49)
- QVAL is wrongly called GQ for indels [\#36](https://github.com/IARCbioinfo/needlestack/issues/36)
- Dockerfile always adds scripts from master branch [\#27](https://github.com/IARCbioinfo/needlestack/issues/27)
- CircleCI deploy.sh doesn't trigger correctly Docker Hub [\#26](https://github.com/IARCbioinfo/needlestack/issues/26)

## [v0.1](https://github.com/IARCbioinfo/needlestack/tree/v0.1) (2015-09-18)
**Implemented enhancements:**

- Choose strand bias filter to apply [\#15](https://github.com/IARCbioinfo/needlestack/issues/15)

**Fixed bugs:**

- Pipeline crashed when no variant is found [\#14](https://github.com/IARCbioinfo/needlestack/issues/14)
- Cutting bed files with zero length regions was not working properly [\#10](https://github.com/IARCbioinfo/needlestack/issues/10)



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
