#!/usr/bin/env bash
#Run as sudo on Apollo2

echo "Downloading latest domains"
wget -O /var/lib/mysql-files/ecod.latest.csv "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt"

sed -i 's/,/;/g' /var/lib/mysql-files/ecod.latest.csv
sed -i 's/\t/,/g' /var/lib/mysql-files/ecod.latest.csv

echo "Uploading data to MYSQL"
mysql -h "130.207.36.76" DESIRE -e "LOAD DATA INFILE '/var/lib/mysql-files/ecod.latest.csv' IGNORE INTO TABLE EcodDomains FIELDS TERMINATED BY ',' ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 5 ROWS;"