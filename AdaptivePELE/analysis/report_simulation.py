from fpdf import FPDF
import re
import sys
import os
import subprocess
from AdaptivePELE.analysis import plotAdaptive
import argparse

ADAPTIVE_PATH = "/run/media/ywest/BAD432E1D432A015/3s3i_out_in_waters/"
IMAGE = "plot_1.png"
CONTACTS = "contacts"
CONTROL_FILE = "originalControlFile_1.conf"
steps_Run, Xcol, Ycol, filename, kind_Print, colModifier, traj_range = 8, 6, 5, "report_", "PRINT_BE_RMSD", 4, None


def arg_parse():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('control_file', type=str, help='adaptive control file')
    args = parser.parse_args()

    return args

def retrieve_metrics(control_file):
    with open(control_file, "r") as f:
        discard = ["random", "constrainAtomToPosition", "sphericalBox"]
        metrics = [line for line in f if (line.strip()).startswith('"type"') and not any(x in line for x in discard)]
    pattern = r'"([A-Za-z0-9_\./\\-]*)"'
    metrics = re.findall(pattern, "".join(metrics))
    metrics = [metric for metric in metrics if metric != "type"]
    return metrics


def write_report(metrics, initial_column=4):

    OUTPUT = 'adaptive_report.pdf'

    pdf = FPDF()
    pdf.add_page() 

    # Title
    pdf.set_font('Arial', 'B', 15)
    pdf.cell(80)
    pdf.cell(30, 10, 'Simulation Report', align='C')

    # Plot Binding SASA
    plot(1+metrics.index("bindingEnergy")+initial_column, 1+metrics.index("sasa")+initial_column, ".", "BE.png")
    pdf.cell(-100)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(10, 49, 'Interaction Energy vs Sasa' +34*"\t" + "Total Energy vs Interaction Energy")
    pdf.image("BE.png", 10, 40, 83)

    # Plot Total Biding
    plot(initial_column,  1+metrics.index("sasa")+initial_column, ".", "total.png")
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image("total.png", 100, 40, 83)

    # Contacts
    create_contact_plot(".")
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image("{}_threshold.png".format(CONTACTS), 10, 110, 83)
    pdf.image("{}_hist.png".format(CONTACTS), 100, 110, 83)
    
    #Plot other metrics agains binding
    doc_position = [-120, 1000]
    for user_metric in ["rmsd", "com_distance", "distanceToPoint"]:
        if user_metric in metrics:
            pdf,doc_position= write_metric(pdf, user_metric, "bindingEnergy", metrics, doc_position)

    #Output report    
    pdf.output(OUTPUT, 'F')


def write_metric(pdf, metric1, metric2, metrics, doc_position, initial_pos=4):
    image_name = "{}.png".format(metric1)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(0)
    pdf.cell(doc_position[0], doc_position[1])
    pdf.cell(8, 10, '{} vs {}'.format(metric1, metric2))
    plot(1+metrics.index(metric1)+initial_pos,  1+metrics.index(metric2)+initial_pos, ".", image_name, zcol=initial_pos + 1+ metrics.index("sasa"))
    pdf.image(image_name, 65, 20, 83)
    return pdf, doc_position
    
def plot(Xcol, Ycol, path, name, zcol=5):
    command = '''gnuplot -e "set terminal png; set output '{}';{}"'''.format(name, plotAdaptive.generatePrintString(8, Xcol, Ycol, "report_", "PRINT_BE_RMSD", zcol, None).strip("\n"))
    os.system(command)
        
def create_contact_plot(path, filenames=CONTACTS):
    command = "python -m AdaptivePELE.analysis.numberOfClusters -f contacts"
    os.system(command)
   
def retrieve_pele_confile(control_file):
    with open(control_file, "r") as f:
        control_file_line = [line for line in f if line.strip().startswith('"controlFile"')][0]
    return control_file_line.split(":")[1].strip().strip('"')

def main(control_file):
    print("Search pele control file")
    pele_conf = retrieve_pele_confile(control_file)
    print("Retrieve metrics")
    metrics = retrieve_metrics(pele_conf)
    print("Build report")
    write_report(metrics)
 

if __name__ == "__main__":
    args = arg_parse()
    main(args.control_file)
