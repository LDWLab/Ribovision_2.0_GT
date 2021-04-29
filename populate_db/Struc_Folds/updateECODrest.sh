#!/usr/bin/env bash

echo "Downloading latest domains"
#wget -O ./ecod.latest.domains.txt "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt"

#sed -i 's/,/;/g' ./ecod.latest.domains.txt
#sed -i 's/\t/,/g' ./ecod.latest.domains.txt

echo "Uploading data to MYSQL"
mysql -h "130.207.36.76" SEREB_test -e "LOAD DATA INFILE './ecod.latest.domains.txt' INTO TABLE SEREB.EcodDomains FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 5 ROWS;"