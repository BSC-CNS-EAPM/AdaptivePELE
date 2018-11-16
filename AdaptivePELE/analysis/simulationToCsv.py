import pandas as pd
import os
import glob

EPOCH = 'Epoch'
TRAJ = 'Traj'


def retrieve_fields(report):
    data = pd.read_csv(report, sep='    ', engine='python')
    return (list(data)[1:2] + list(data)[3:])

def gather_reports():
    reports = glob.glob(os.path.join("*/*report*"))
    reports = glob.glob(os.path.join("*report*")) if not reports else reports
    reports = filter_non_numerical_folders(reports)
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")
    return reports

def filter_non_numerical_folders(reports, numfolders=True):
    """
    Filter non numerical folders among
    the folders to parse
    """
    if(numfolders):
        new_reports = [report for report in reports if(os.path.basename(os.path.dirname(report)).isdigit())]
        return new_reports
    else:
        return reports


def fill_data(reports, df):

    for report in reports:
       print(report)
       data = pd.read_csv(report, sep='    ', engine='python')
       data = data.drop(data.columns[0], axis=1)
       data = data.drop(data.columns[1], axis=1)
       size = data.shape[0]
       data.insert(0, EPOCH, [os.path.dirname(report)]*size)
       data.insert(1, TRAJ, [os.path.basename(report).split("_")[-1]]*size)
       df = pd.concat([df, data], axis=0)
    return df
       

if __name__ == "__main__":
   epochs = [folder for folder in glob.glob("./*/") if folder.isdigit()]
   reports = gather_reports()
   fields = retrieve_fields(reports[0])
   df = pd.DataFrame({EPOCH : [], TRAJ : [] })
   for field in fields:
       df[field] = []
   df = fill_data(reports, df)
   df.to_csv("report_summary.csv", index=False, decimal=".", float_format='%.2f')
