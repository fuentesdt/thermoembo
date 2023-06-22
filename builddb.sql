-- cat builddb.sql  | sqlite3 database/ctkDICOM.sql

.output database/datakey.csv 
select "" ImageType, pt.PatientID,sd.StudyInstanceUID,sd.StudyDate,se.SeriesInstanceUID,se.SeriesDescription,pt.PatientID || "/" ||sd.StudyInstanceUID || "/" || sd.StudyDate || "/" || se.SeriesInstanceUID as filepath
from patients pt
join studies sd on pt.UID=sd.PatientsUID
join series se  on se.StudyInstanceUID=sd.StudyInstanceUID;
