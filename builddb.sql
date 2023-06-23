-- cat builddb.sql  | sqlite3 database/ctkDICOM.sql
-- sqlite3 -init builddb.sql  database/ctkDICOM.sql

attach database ':memory:' as tmp;
attach database 'database/ctkDICOMTagCache.sql' as tag;
create table tmp.flagdata  as
select CASE WHEN se.seriesDescription like "%late%art%3%" THEN 'cta'
     WHEN se.SeriesDescription like "%venous%3%"  THEN 'triplephase'
     ELSE NULL END AS ImageType, pt.PatientID,sd.StudyInstanceUID,sd.StudyDate,se.SeriesInstanceUID,se.SeriesDescription,im.SOPInstanceUID,pt.PatientID || "/" ||sd.StudyInstanceUID || "/" || sd.StudyDate || "/" || se.SeriesInstanceUID as filepath
from patients pt
join studies sd on pt.UID=sd.PatientsUID
join series se  on se.StudyInstanceUID=sd.StudyInstanceUID 
join images im  on im.SeriesInstanceUID=se.SeriesInstanceUID
group by se.SeriesInstanceUID;

-- slice thickness not in tags :(
--select fd.*,tc.* 
--from tmp.flagdata fd
--join tag.TagCache tc on fd.SOPInstanceUID = tc.SOPInstanceUID 
--where tc.Tag = '0018,0050' ;
--
--select * from tag.TagCache tc where tc.Tag = '0018,0050' ;

.headers on
.output database/datakey.csv 
select * from tmp.flagdata;
