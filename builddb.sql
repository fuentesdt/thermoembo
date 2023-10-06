-- cat builddb.sql  | sqlite3 database/ctkDICOM.sql
-- sqlite3 -init builddb.sql  database/ctkDICOM.sql

attach database ':memory:' as tmp;
attach database 'database/ctkDICOMTagCache.sql' as tag;
create table tmp.flagdata  as
select CASE WHEN se.seriesDescription like "%arterial%3%I%" THEN 'art'
     WHEN se.seriesDescription like "%late%art%3%I%" THEN 'cta'
     WHEN se.SeriesDescription like "%venous%3%I%"    THEN 'ven'
     WHEN se.SeriesDescription like "%non%contrast%3%I%"  THEN 'pre'
     ELSE NULL END AS ImageType, pt.PatientID,sd.StudyInstanceUID,sd.StudyDate,se.SeriesInstanceUID,se.SeriesDescription,se.SeriesNumber,im.SOPInstanceUID,pt.PatientID || "/" ||sd.StudyDate || "/" ||sd.StudyInstanceUID || "/" ||  se.SeriesInstanceUID as Filename
from patients pt
join studies sd on pt.UID=sd.PatientsUID
join series se  on se.StudyInstanceUID=sd.StudyInstanceUID 
join images im  on im.SeriesInstanceUID=se.SeriesInstanceUID
group by se.SeriesInstanceUID;



create table tmp.treatmentdates(mrn text primary key, tumorinduction date, treatment date, necropsy date  );
INSERT INTO tmp.treatmentdates (mrn,tumorinduction,treatment,necropsy)
VALUES 
('ZPAF23S018','2023-05-17','2023-05-31','2023-06-07'),
('ZPAF23S019','2023-03-29','2023-04-12','2023-04-18'),
('ZPAF23S020','2023-04-05','2023-04-18','2023-04-26'),
('ZPAF23S021','2023-04-19','2023-05-13','2023-05-10'),
('ZPAF23S022','2023-05-03','2023-05-17','2023-05-17'),
('ZPAF23S023','2023-05-10','2023-05-24','2023-05-31');


create table tmp.widestudy  as
select fg.PatientID,
       fg.StudyDate,
            group_concat(CASE WHEN ImageType = 'pre' THEN fg.SeriesDescription       ELSE NULL END)  PREDescription,
            group_concat(CASE WHEN ImageType = 'art' THEN fg.SeriesDescription       ELSE NULL END)  ARTDescription,
            group_concat(CASE WHEN ImageType = 'ven' THEN fg.SeriesDescription       ELSE NULL END)  VENDescription,
            group_concat(CASE WHEN ImageType = 'cta' THEN fg.SeriesDescription       ELSE NULL END)  CTADescription,
       fg.StudyInstanceUID,
            group_concat(CASE WHEN ImageType = 'pre' THEN fg.SeriesInstanceUID       ELSE NULL END)  PREUID,
            group_concat(CASE WHEN ImageType = 'pre' THEN fg.Filename                ELSE NULL END)  PREFilename,
            group_concat(CASE WHEN ImageType = 'art' THEN fg.SeriesInstanceUID       ELSE NULL END)  ARTUID,
            group_concat(CASE WHEN ImageType = 'art' THEN fg.Filename                ELSE NULL END)  ARTFilename,
            group_concat(CASE WHEN ImageType = 'ven' THEN fg.SeriesInstanceUID       ELSE NULL END)  VENUID,
            group_concat(CASE WHEN ImageType = 'ven' THEN fg.Filename                ELSE NULL END)  VENFilename,
            group_concat(CASE WHEN ImageType = 'cta' THEN fg.SeriesInstanceUID       ELSE NULL END)  CTAUID,
            group_concat(CASE WHEN ImageType = 'cta' THEN fg.Filename                ELSE NULL END)  CTAFilename
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
select td.tumorinduction,td.treatment,td.necropsy,ws.* from tmp.widestudy ws left join tmp.treatmentdates td on ws.patientid=td.mrn order by ws.PatientID,ws.StudyDate;
.quit
