-- MySQL code to create a database of cells, cell information from Darmanis et al. and SC3 assigned clusters (k=9) from PigeonPigs
-- requires the count_for_MySQL.csv produced by count_table_maker.py, the downloaded SRA Run Table from Darmanis with irrelevant columns removed and the SC3 k9 clusters per cell from PigeonPigs
--run this is mysql on the server to produce a VIEW that users can query

DROP DATABASE if exists dbGroupA;
CREATE DATABASE dbGroupA;
USE dbGroupA;

CREATE TABLE RunTable (
    Run VARCHAR(10),
    SRA_Sample VARCHAR(10),
    Sample_Name VARCHAR(10),
    age VARCHAR(20),
    c1_chip_id VARCHAR(20),
    cell_type VARCHAR(20),
    experiment_sample_name VARCHAR(10),
    tissue VARCHAR(20),
    PRIMARY KEY (Run)
);

CREATE TABLE k_cluster (
    Run VARCHAR(10),
    k9_cluster INT,
    PRIMARY KEY (Run)
);

CREATE TABLE counts (
    Run VARCHAR(10),
    Gene VARCHAR(40),
	Gene_Count INT NOT NULL,
    PRIMARY KEY (Run, Gene)
);

LOAD DATA LOCAL INFILE '/home/graham/Downloads/Wobsite/MySQL/Table_for_SQL.csv' INTO TABLE RunTable
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;

LOAD DATA LOCAL INFILE '/home/graham/Downloads/Wobsite/MySQL/count_for_MySQL.csv' INTO TABLE counts
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

LOAD DATA LOCAL INFILE '/home/graham/Downloads/Wobsite/MySQL/k9_clusters.csv' INTO TABLE k_cluster
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

CREATE VIEW website AS SELECT * FROM counts NATURAL JOIN RunTable NATURAL JOIN k_cluster ORDER BY Gene;
