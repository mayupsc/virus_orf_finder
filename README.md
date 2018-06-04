# virus_orf_finder
A web app based on shiny, DECIPHER, NCBI ORFfinder and BLAST
# Virus ORF Finder workflow

### 1. preparatory work

pull docker image
```
docker pull yumapsc/virus_orf_finder

```
clone the repository

```
git clone https://github.com/mayupsc/virus_orf_finder
```

run docker image

```
docker run --name shiny -d -p 3838:3838 -v viral_genome:/srv/shiny-server/ yumapsc/virus_orf_finder
```

### 2. Upload Virus genome files
###

### 3. Set Minimal ORF length for ORF finder (default:50)
###

### 4. Choose interest virus(default:AJ489258.1) and set Expectation value(E) threshold for BLAST(default:1e-10)
###

### 5. Click 'Calculator' menu and wait for two minutes
###

### 6. Click 'One ORF' submenu and then you will see hits alignment and distribution.You can select the ORF that you are 
###
###    interested in.
