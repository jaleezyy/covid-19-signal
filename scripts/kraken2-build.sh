## download taxonomy
kraken2-build --download-taxonomy --db db --threads 20

## download db
kraken2-build --download-library viral --db db --threads 20
#kraken2-build --download-library human --db db --threads 20

## clean before building db
rm db/hash.k2d  db/opts.k2d  db/seqid2taxid.map  db/taxo.k2d

## build
kraken2-build --build --threads 20 --db db

## clean intermediate files in db
kraken2-build --clean --threads 20 --db db
