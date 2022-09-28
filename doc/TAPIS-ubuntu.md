# Long-reads processing with TAPIS on Ubuntu 20.04 docker container


# Docker container preparation


```
ls /Users/jsun/Desktop/silk
## input.fa chromosome.fa

docker pull ubuntu:20.04
docker run -v /Users/jsun/Desktop/silk:/home/ubuntu -it ubuntu:20.04 bash
```


# Environment preparation

```
pwd
# /home/ubuntu

# install software onto Ubuntu
apt update
apt install vim curl wget build-essential automake cmake \
            libpng-dev libfreetype6-dev libfontconfig1-dev xclip \
            liblzma-dev libncurses5-dev libbz2-dev libcurl4-nss-dev \
            git python python-dev 
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
python get-pip.py


# install samtools
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -jxvf htslib-1.16.tar.bz2
tar -jxvf samtools-1.16.1.tar.bz2
cd htslib-1.16
./configure
make
make install
cd ..
cd samtools-1.16.1
./configure
make
make install
cd ..


# install GMAP v2016-06-30
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2016-06-30.v5.tar.gz
tar xzvf gmap-gsnap-2016-06-30.v5.tar.gz
cd gmap-2016-06-30
./configure
make
make install
cd ..


# install SpliceGrapher v0.2.5
# manually download SpliceGrapher-0.2.5.tgz and share the file to the Docker container
tar xzvf SpliceGrapher-0.2.5.tgz
cd SpliceGrapher-0.2.5
python setup.py install
cd ..


# install TAPIS
git clone git@github.com:jsun/TAPIS.git
cd TAPIS
pip install -r requirements.txt
python setup.py build
python setup.py install
cd ..
```

# Data analysis

```
pwd
# /home/ubuntu

mkdir gmap_index
gmap_build -D gmap_index -d chromosome chromosome.fa
alignPacBio.py -v -o tapis_outputs gmap_index chromosome chromosome.fa input.fa
```



