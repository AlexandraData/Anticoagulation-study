# Load required libraries
library(DatabaseConnector)
library(Eunomia)
library(dplyr)

# 1. Create a permanent SQLite database file inside your project's data folder
db_path <- file.path("data", "eunomia_cdm.sqlite")

# getDatabaseFile automatically extracts the downloaded GiBleed zip into the SQLite file
getDatabaseFile(
  datasetName = "GiBleed", 
  cdmVersion = "5.3",
  databaseFile = db_path,
  overwrite = TRUE
)

# 2. Configure the database connection using DatabaseConnector
connectionDetails <- createConnectionDetails(
  dbms = "sqlite",
  server = db_path
)

conn <- connect(connectionDetails)

# 3. Verify the connection and check patient count
patient_count <- querySql(conn, "SELECT COUNT(*) AS total_patients FROM person;")
print(patient_count)

# 4. View top 10 prescribed drugs in the dataset
top_drugs_sql <- "
SELECT c.concept_name, c.concept_id, COUNT(*) AS exposure_count
FROM drug_exposure de
JOIN concept c ON de.drug_concept_id = c.concept_id
GROUP BY c.concept_name, c.concept_id
ORDER BY exposure_count DESC
LIMIT 10;
"
top_drugs <- querySql(conn, top_drugs_sql)
print(top_drugs)

#---------------------------------------------------------------
# Query anticoagulant & antiplatelet drug concepts in the dataset
anticoag_check_sql <- "
SELECT c.concept_name, c.concept_id, COUNT(*) AS patient_count
FROM drug_exposure de
JOIN concept c ON de.drug_concept_id = c.concept_id
WHERE LOWER(c.concept_name) LIKE '%warfarin%'
   OR LOWER(c.concept_name) LIKE '%apixaban%'
   OR LOWER(c.concept_name) LIKE '%rivaroxaban%'
   OR LOWER(c.concept_name) LIKE '%aspirin%'
   OR LOWER(c.concept_name) LIKE '%clopidogrel%'
   OR LOWER(c.concept_name) LIKE '%celecoxib%'
   OR LOWER(c.concept_name) LIKE '%diclofenac%'
GROUP BY c.concept_name, c.concept_id
ORDER BY patient_count DESC;
"

anticoag_summary <- querySql(conn, anticoag_check_sql)
print(anticoag_summary)