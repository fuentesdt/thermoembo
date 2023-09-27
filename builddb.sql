-- cat builddb.sql  | sqlite3 database/ctkDICOM.sql
-- sqlite3 -init builddb.sql  database/ctkDICOM.sql

attach database ':memory:' as tmp;
attach database 'database/ctkDICOMTagCache.sql' as tag;
create table tmp.flagdata  as
select CASE WHEN se.seriesDescription like "%arterial%3%I%40%" THEN 'art'
     WHEN se.seriesDescription like "%late%art%3%I%40%" THEN 'cta'
     WHEN se.SeriesDescription like "%venous%3%I%"    THEN 'ven'
     WHEN se.SeriesDescription like "%non%contrast%"  THEN 'pre'
     ELSE NULL END AS ImageType, pt.PatientID,sd.StudyInstanceUID,sd.StudyDate,se.SeriesInstanceUID,se.SeriesDescription,im.SOPInstanceUID,pt.PatientID || "/" ||sd.StudyDate || "/" ||sd.StudyInstanceUID || "/" ||  se.SeriesInstanceUID as Filename
from patients pt
join studies sd on pt.UID=sd.PatientsUID
join series se  on se.StudyInstanceUID=sd.StudyInstanceUID 
join images im  on im.SeriesInstanceUID=se.SeriesInstanceUID
group by se.SeriesInstanceUID;

create table tmp.widestudy  as
select fg.PatientID,fg.StudyInstanceUID, fg.StudyDate,
            max(CASE WHEN ImageType = 'pre' THEN fg.SeriesDescription       ELSE NULL END)  PREDescription,
            max(CASE WHEN ImageType = 'pre' THEN fg.SeriesInstanceUID       ELSE NULL END)  PREUID,
            max(CASE WHEN ImageType = 'pre' THEN fg.Filename                ELSE NULL END)  PREFilename,
            max(CASE WHEN ImageType = 'art' THEN fg.SeriesDescription       ELSE NULL END)  ARTDescription,
            max(CASE WHEN ImageType = 'art' THEN fg.SeriesInstanceUID       ELSE NULL END)  ARTUID,
            max(CASE WHEN ImageType = 'art' THEN fg.Filename                ELSE NULL END)  ARTFilename,
            max(CASE WHEN ImageType = 'ven' THEN fg.SeriesDescription       ELSE NULL END)  VENDescription,
            max(CASE WHEN ImageType = 'ven' THEN fg.SeriesInstanceUID       ELSE NULL END)  VENUID,
            max(CASE WHEN ImageType = 'ven' THEN fg.Filename                ELSE NULL END)  VENFilename,
            max(CASE WHEN ImageType = 'cta' THEN fg.SeriesDescription       ELSE NULL END)  CTADescription,
            max(CASE WHEN ImageType = 'cta' THEN fg.SeriesInstanceUID       ELSE NULL END)  CTAUID,
            max(CASE WHEN ImageType = 'cta' THEN fg.Filename                ELSE NULL END)  CTAFilename
from tmp.flagdata  fg
where fg.ImageType is not null 
GROUP BY    fg.StudyInstanceUID;
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
.output database/datawide.csv 
select * from tmp.widestudy;
.quit
