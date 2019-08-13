'''
----------------------------------------------------------------------------------
|             application2.py - Data visualization Flask application             |
----------------------------------------------------------------------------------
| Kristoph Naggert | Mount Desert Island Biological Laboratory | Oberlin College |
----------------------------------------------------------------------------------
|                           Friday, June 18th, 2019                              |
----------------------------------------------------------------------------------
'''

'''
Import required libraries and specific functions.
'''

from flask import Flask, render_template, request
from wtforms import Form, TextField, validators
from wtforms.validators import DataRequired, Optional
from datetime import datetime
import sys, time, requests, subprocess, os, re
import os.path
from os import path
import pandas as pd
import numpy as np
from math import pi

'''
Import Bokeh Functions
'''
from bokeh.io import show
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import LinearColorMapper, BasicTicker, PrintfTickFormatter, ColorBar, ContinuousTicker, CheckboxButtonGroup, CheckboxGroup, CustomJS
from bokeh.models import Legend, LegendItem, Range1d, HoverTool, RedoTool, UndoTool
from bokeh.models.widgets import Tabs, Panel
from bokeh.embed import components
from bokeh.layouts import layout, widgetbox, column, row, gridplot
from bokeh.palettes import Viridis5, Magma256, PiYG
from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead


application = Flask(__name__)


'''
Initialize Global Varaibles
'''

gene_from_genome = ''
biological_identifier = ''
upstreamBuf = ''
downstreamBuf = ''
user_input_seq = ''
seq = ''
file_name = ''
script = ''
div = ''
script2 = ''
div2 = ''



class RequiredIf(DataRequired):

    '''
    Validator which makes a field required if another field is set and has a truthy value.

    Sources:
        - http://wtforms.simplecodes.com/docs/1.0.1/validators.html
        - http://stackoverflow.com/questions/8463209/how-to-make-a-field-conditionally-optional-in-wtforms
        - https://gist.github.com/devxoul/7638142#file-wtf_required_if-py
    '''

    field_flags = ('requiredif',)

    def __init__(self, message=None, *args, **kwargs):
        super(RequiredIf).__init__()
        self.message = message
        self.conditions = kwargs

    # field is requiring that name field in the form is data value in the form
    def __call__(self, form, field):
        for name, data in self.conditions.items():
            other_field = form[name]
            if other_field is None:
                raise Exception('no field named "%s" in form' % name)
            if other_field.data == data and not field.data:
                DataRequired.__call__(self, form, field)
            Optional()(form, field)



class InputForm(Form):
    gene = TextField(label = 'Enter Standard Yeast Gene ID or Gene Name', validators = [validators.optional()])
    sequence = TextField(label = 'Enter User Generated Sequence', validators = [RequiredIf(gene='')])
    identifier = TextField(label = 'Enter a Unique Gene Identifier for User Inputed Sequence, Ex. >YLR115W. chromosome:R64-1-1:XII:377485:380664:1', validators = [RequiredIf(gene = '')])
    upstream_buffer = TextField(label = 'Enter Length for 5-Prime End Buffer', validators = [RequiredIf(gene='')], default='100')
    downstream_buffer = TextField(label = 'Enter Length for 3-Prime End Buffer', validators = [RequiredIf(gene='')], default='500')



@application.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():

        global gene_from_genome
        global biological_identifier
        global upstreamBuf
        global downstreamBuf
        global user_input_seq

        gene_from_genome = form.gene.data
        upstreamBuf = form.upstream_buffer.data
        downstreamBuf = form.downstream_buffer.data
        user_input_seq = form.sequence.data
        biological_identifier = form.identifier.data


        now = datetime.now()
        current_timestamp = str(datetime.timestamp(now))
        print("timestamp =", current_timestamp)

        global file_name

        if (len(gene_from_genome) != 0 and len(user_input_seq) == 0):
            gene_from_genome = convert_to_symbol(gene_from_genome)
            print(gene_from_genome)
            if (gene_from_genome == 'dne'):
                print('heyo, made it here')
                return render_template('error2.html')
            else:
                seq_check = False
                file_name = gene_from_genome+'_'+upstreamBuf+'_'+downstreamBuf+'_'+current_timestamp+'.fa'
                fasta_file = open('/var/kristoph_flask/outfiles/'+file_name, 'w')
                fasta_file.write(get_seq(gene_from_genome, upstreamBuf, downstreamBuf))
                fasta_file.close()
        else:
            seq_check = True
            file_name = 'userseq'+'_'+current_timestamp+'.fa'
            fasta_file = open('/var/kristoph_flask/outfiles/'+file_name, 'w')
            fasta_file.write(biological_identifier+'\n')
            fasta_file.write(user_input_seq)
            fasta_file.close()

        create_batchScript(file_name)
        
        processed_file = file_name+'.pos.txt'

        print(processed_file)
        path_to_file = '/var/www/vhosts/knaggert-vm.mdibl.net/flask_project/paHMM/testHMM/Export/'+processed_file

        while not path.exists(path_to_file):
            time.sleep(1)

        
        raw_data = pd.read_csv(path_to_file, sep = '\s+', skiprows = 1, header = None, names = ['Base', 'Position', 'e1', 'e2', 'e3', 'pASite', 'e4'])

        while(raw_data.empty):
            print(raw_data.empty)
            time.sleep(1)
            raw_data = pd.read_csv(path_to_file, sep = '\s+', skiprows = 1, header = None, names = ['Base', 'Position', 'e1', 'e2', 'e3', 'pASite', 'e4'])

        TOOLS = "hover, save, box_zoom, pan, undo, redo, reset, wheel_zoom, tap"

        print("Check 1")
        col1, col2, col3 = gridded_plots(raw_data, TOOLS, gene_from_genome, upstreamBuf)
        
        print("Check 2")
        heatmap = Intensity_Plot(raw_data, TOOLS)

        print("Check 3")
        l1 = gridplot([[col3]])
        l2 = gridplot([[heatmap]])

        print("Check 4")
        tab1 = Panel(child=l1, title="Line")
        tab2 = Panel(child=l2, title="heatmap")
        
        print("Check 5")
        tab = Tabs(tabs=[ tab1, tab2 ])

        print("Check 6")
        l3 = gridplot([[col1, tab]])
        l4 = gridplot([[col2, tab]])

        print("Check 7")
        t1 = Panel(child=l3, title="Independent")
        t2 = Panel(child=l4, title="Staged")

        print("Check 8")
        tabs = Tabs(tabs=[ t1, t2 ])

        global script
        global div

        script, div = components(tabs)
        
        print('complete, no errors')

    return render_template('query2.html', form = form, gene_from_genome = gene_from_genome, upstreamBuf = upstreamBuf, downstreamBuf = downstreamBuf, user_input_seq = user_input_seq, seq = seq, script = script, div = div)


application.route('/about')
def about():
    return render_template('about.html')

application.route('/datadownload')
def downloaddata():
    return render_template('download.html')

application.route('/troubleshooting')
def troubleshooting():
    return render_template('troubleshooting.html')

application.route('/literature')
def lit():
    return render_template('literature.html')


def create_batchScript(file_name):

    trunc_file = file_name[0:-3]
    print(trunc_file)
    test_file_name = 'testScript_'+trunc_file+'.txt'
    print(test_file_name)
    batchScript = open('/var/kristoph_flask/outfiles/testScripts/'+test_file_name, 'w')
    batchScript.write('/var/www/vhosts/knaggert-vm.mdibl.net/flask_project/paHMM/testHMM\n')
    batchScript.write('load\n')
    batchScript.write('yeastHMM\n')
    batchScript.write('apply\n')
    batchScript.write('yeastP\n')
    batchScript.write('3 22 41 55 64 -1\n')
    batchScript.write('/var/kristoph_flask/outfiles/fa_files/'+file_name+'\n')
    batchScript.write('exit')
    batchScript.close()

    return

def get_seq(gene_from_genome, upstreamBuf, downstreamBuf):

    '''
    Sourced from: https://rest.ensembl.org/documentation/info/sequence_id
    '''

    server = "https://rest.ensembl.org"

    ext = '/sequence/id/'+gene_from_genome+'?expand_5prime='+upstreamBuf+';expand_3prime='+downstreamBuf

    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    seq = r.text

    return seq



def convert_to_symbol(gene_from_genome):

    ''' Read in the data '''

    df = pd.read_csv('/var/kristoph_flask/data/results.tsv', delimiter= '\s+', header = None, names = ['Sys_Name','Std_Name'])

    gene_from_genome = gene_from_genome.upper()

    if (any(df['Sys_Name'] == gene_from_genome)):
        gene_from_genome = gene_from_genome
        return gene_from_genome
    elif ((all(df['Sys_Name'] != gene_from_genome) and (all(df['Std_Name'] != gene_from_genome)))):
        gene_from_genome = 'dne'
        return gene_from_genome
    else:
        index = next(iter(df[df.Std_Name==gene_from_genome].index), 'No match.')
        gene_from_genome = df.Sys_Name[index]
        return gene_from_genome

    '''
    Currently has no way to catch and inform the user of the application
    if the standard gene name or systematic gene name do not exist in the
    local database.

    Also need to create a catch if there is no input, becuase there are certain
    entries in the gene list (results.tsv) whose Std_Name entry is ""and the
    application might break if the program tries to find the gene name "" when a
    user sequence is input.
    '''

def Intensity_Plot(data, TOOLS):

    dataI = data
    
    columns = ['e1', 'e2', 'e3', 'pASite', 'e4']

    array_1 = []

    for k in range(len(dataI)):
        if (k+6 >= len(dataI)):
            array_1.append("N/A")
        else:
            s = ""
            list1 = []
            for l in range(k, k+6):
                list1.append(dataI.Base[l])
            s = s.join(list1)
            array_1.append(s)

    sequences = []
    for n in range(len(dataI)*len(columns)):
        sequences.append(n)

    temporary = 0
    for q in range(len(array_1)):
        for d in range(len(columns)):
            sequences[temporary] = array_1[q]
            temporary = temporary+1

    Sequences = pd.DataFrame(sequences, columns = ['Sequence'])

    dataI['Position'] = dataI['Position'].astype(str)
    dataI = dataI.set_index('Position')

    dataI = dataI.drop(columns = ['Base'])

    dataI.columns.name = 'Values'
    Position = list(dataI.index)
    Position = list(map(int, Position))
    Values = list(dataI.columns)

    df_I = pd.DataFrame(dataI.stack(), columns = ['score']).reset_index()

    df_complete = pd.concat([df_I, Sequences], axis = 1)

    print(df_complete.head(20))
    
    Magma256.reverse()
    mapper = LinearColorMapper(palette = Magma256, high = df_complete.score.max(), low = df_complete.score.min())

    p = figure(title = "Heatmap of Sites", x_range = (0,len(Position)), y_range = Values, 
               x_axis_location = "above", sizing_mode = 'stretch_both',
               tools = TOOLS, toolbar_location = 'below',
               tooltips = [('Position','@Position'), ('Score', '@score'), ('Sequence', '@Sequence')])



    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "10pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3

    p.rect(x = 'Position', y = 'Values', width = 1, height = 1,
           source = df_complete,
           fill_color = {'field':'score', 'transform' : mapper},
           line_color = None)

    color_bar = ColorBar(color_mapper = mapper, major_label_text_font_size = "10pt",
                         ticker = BasicTicker(),
                         label_standoff = 6, border_line_color = None, location = (0, 0))

    p.add_layout(color_bar, 'right')

    return(p)


def gridded_plots(raw_data, TOOLS, gene_from_genome, upstreamBuf):

    dataS = raw_data

    '''
    Stacked Plots (Left-Hand Side)
    '''

    array = []

    for k in range(len(dataS)):
        if (k+6 >= len(dataS)):
            array.append("N/A")
        else:
            s = ""
            list_1 = []
            for l in range(k, k+6):
                list_1.append(dataS.Base[l])
            s = s.join(list_1)
            array.append(s)

    Seq = pd.DataFrame(array, columns = ['Seq'])
    df_comp = pd.concat([dataS, Seq], axis = 1)

    a = 1.5
    b = 18.5

    temp = []
    for q in range(len(df_comp)):
        temp.append(1 - 1/(1 + 2**(a*(df_comp.pASite[q] - b))))

    df_comp['pASiteT'] = temp

    temp2 = []

    for s in range(len(df_comp)):
        if (s == 0):
            temp2.append(df_comp['pASiteT'][0])
        else:
            temp2.append(df_comp['pASiteT'][s]+temp2[s-1])

    df_comp['cumulative'] = temp2

    last = len(df_comp)
    L_entry = df_comp.cumulative[last-1]

    df_comp['cumulative'] = df_comp['cumulative']/L_entry
    df_comp['pASiteT'] = df_comp['pASiteT']/L_entry

    source = ColumnDataSource(df_comp)

    ''' Individual Line Graph for  e1 '''

    s1 = figure(title = 'e1', toolbar_location = "right", tools = TOOLS, tooltips = [('Position', '@Position'), ('Score', '@e1'), ('Sequence', '@Seq')])

    s1.line('Position', 'e1', line_width = 2, color = Viridis5[0], source = source)

    s1.xaxis.visible = False

    ''' Individual Line Graph for  e2 '''

    s2 = figure(title = 'e2', x_range = s1.x_range, y_range = s1.y_range, tools = TOOLS, tooltips = [('Position', '@Position'), ('Score', '@e2'), ('Sequence', '@Seq')])

    s2.line('Position', 'e2', line_width = 2, color = Viridis5[1], source = source)

    s2.xaxis.visible = False

    ''' Individual Line Graph for  e3 '''

    s3 = figure(title = 'e3',x_range = s1.x_range, y_range = s1.y_range, tools = TOOLS, tooltips = [('Position', '@Position'), ('Score', '@e3'), ('Sequence', '@Seq')])

    s3.line('Position', 'e3', line_width = 2, color = Viridis5[2], source = source)

    s3.xaxis.visible = False

    ''' Individual Line Graph for  pA Site '''

    s4 = figure(title = 'pA Site', x_range = s1.x_range, y_range = s1.y_range, tools = TOOLS,
                tooltips = [('Position', '@Position'), ('Score', '@pASite'), ('Sequence', '@Seq')])

    s4.line('Position', 'pASite', line_width = 2, color = Viridis5[3], source = source)

    s4.xaxis.visible = False

    ''' Individual Line Graph for  e4 '''

    s5 = figure(title='e4',x_range = s1.x_range, y_range = s1.y_range, tools = TOOLS, tooltips = [('Position', '@Position'), ('Score', '@e4'), ('Sequence', '@Seq')])

    s5.line('Position', 'e4', line_width = 2, color = Viridis5[4], source = source)

    s5.yaxis.axis_label = "Score"
    s5.xaxis.axis_label = "Position"

    ''' PCF11 '''

    file = open('/var/kristoph_flask/data/pcf11_cumPa.txt', 'r')

    for line in file:
        if re.search(gene_from_genome, line):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneData.txt', 'a+') as f:
                f.write(line)

    file2 = open('/var/kristoph_flask/data/pcf11_paProb.txt', 'r')

    for line2 in file2:
        if re.search(gene_from_genome, line2):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneProb.txt', 'a+') as fa:
                fa.write(line2)


    data = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneData.txt', sep = '\s+', header = None, names = ['#gene', 'chromo', 'position', 'strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount','AGcount',
                                                                                                                                    'Acount', 'Gcount', 'DHch01_20nt_Ttrim', 'DHch02_20nt_Ttrim', 'DHch03_20nt_Ttrim', 'DHch04_20nt_Ttrim', 'DHch05_20nt_Ttrim', 'DHch06_20nt_Ttrim', 'DHch07_20nt_Ttrim', 'DHch08_20nt_Ttrim', 'all'])

    df = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneProb.txt', sep = '\s+', header = None, names = ['#gene', 'chromo', 'position', 'strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount','AGcount',
                                                                                                                                  'Acount', 'Gcount', 'DHch01_20nt_Ttrim', 'DHch02_20nt_Ttrim', 'DHch03_20nt_Ttrim', 'DHch04_20nt_Ttrim', 'DHch05_20nt_Ttrim', 'DHch06_20nt_Ttrim', 'DHch07_20nt_Ttrim', 'DHch08_20nt_Ttrim', 'all'])


    data = data.drop(['strand', 'avgGcount','AGcount', 'Acount', 'Gcount', 'DHch01_20nt_Ttrim', 'DHch02_20nt_Ttrim', 'DHch03_20nt_Ttrim', 'DHch04_20nt_Ttrim', 'DHch05_20nt_Ttrim', 'DHch06_20nt_Ttrim', 'DHch07_20nt_Ttrim', 'DHch08_20nt_Ttrim'], axis = 1)
    df = df.drop(['strand', 'avgGcount','AGcount', 'Acount', 'Gcount', 'DHch01_20nt_Ttrim', 'DHch02_20nt_Ttrim', 'DHch03_20nt_Ttrim', 'DHch04_20nt_Ttrim', 'DHch05_20nt_Ttrim', 'DHch06_20nt_Ttrim', 'DHch07_20nt_Ttrim', 'DHch08_20nt_Ttrim'], axis = 1)


    ''' decay1 '''


    file3 = open('/var/kristoph_flask/data/decay_cumPa.txt', 'r')

    for line3 in file3:
        if re.search(gene_from_genome, line3):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneData.txt', 'a+') as f:
                f.write(line3)

    file4 = open('/var/kristoph_flask/data/decay_paProb.txt', 'r')

    for line4 in file4:
        if re.search(gene_from_genome, line4):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneProb.txt', 'a+') as fa:
                fa.write(line4)

    data2 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneData.txt', sep = '\s+', header = None, names = ['#gene', 'chromo', 'position', 'strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e1_filtered_uniq30nt', 'UA05e1_filtered_uniq30nt', 'UA10e1_filtered_uniq30nt', 'UA20e1_filtered_uniq30nt', 'UA40e1_filtered_uniq30nt', 'UB00e1_filtered_uniq30nt', 'UB05e1_filtered_uniq30nt', 'UB10e1_filtered_uniq30nt', 'UB20e1_filtered_uniq30nt','UB40e1_filtered_uniq30nt','YA00e1_filtered_uniq30nt','YA05e1_filtered_uniq30nt','YA10e1_filtered_uniq30nt','YA20e1_filtered_uniq30nt','YA40e1_filtered_uniq30nt','YB00e1_filtered_uniq30nt','YB05e1_filtered_uniq30nt','YB10e1_filtered_uniq30nt','YB20e1_filtered_uniq30nt','YB40e1_filtered_uniq30nt','all'])

    df2 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneProb.txt', sep = '\s+', header = None, names = ['#gene', 'chromo', 'position', 'strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e1_filtered_uniq30nt', 'UA05e1_filtered_uniq30nt', 'UA10e1_filtered_uniq30nt', 'UA20e1_filtered_uniq30nt', 'UA40e1_filtered_uniq30nt', 'UB00e1_filtered_uniq30nt', 'UB05e1_filtered_uniq30nt', 'UB10e1_filtered_uniq30nt', 'UB20e1_filtered_uniq30nt','UB40e1_filtered_uniq30nt','YA00e1_filtered_uniq30nt','YA05e1_filtered_uniq30nt','YA10e1_filtered_uniq30nt','YA20e1_filtered_uniq30nt','YA40e1_filtered_uniq30nt','YB00e1_filtered_uniq30nt','YB05e1_filtered_uniq30nt','YB10e1_filtered_uniq30nt','YB20e1_filtered_uniq30nt','YB40e1_filtered_uniq30nt','all'])

    data2 = data2.drop(['strand', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e1_filtered_uniq30nt', 'UA05e1_filtered_uniq30nt', 'UA10e1_filtered_uniq30nt', 'UA20e1_filtered_uniq30nt', 'UA40e1_filtered_uniq30nt', 'UB00e1_filtered_uniq30nt', 'UB05e1_filtered_uniq30nt', 'UB10e1_filtered_uniq30nt', 'UB20e1_filtered_uniq30nt','UB40e1_filtered_uniq30nt','YA00e1_filtered_uniq30nt','YA05e1_filtered_uniq30nt','YA10e1_filtered_uniq30nt','YA20e1_filtered_uniq30nt','YA40e1_filtered_uniq30nt','YB00e1_filtered_uniq30nt','YB05e1_filtered_uniq30nt','YB10e1_filtered_uniq30nt','YB20e1_filtered_uniq30nt','YB40e1_filtered_uniq30nt'], axis=1)

    df2 =  df2.drop(['strand', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e1_filtered_uniq30nt', 'UA05e1_filtered_uniq30nt', 'UA10e1_filtered_uniq30nt', 'UA20e1_filtered_uniq30nt', 'UA40e1_filtered_uniq30nt', 'UB00e1_filtered_uniq30nt', 'UB05e1_filtered_uniq30nt', 'UB10e1_filtered_uniq30nt', 'UB20e1_filtered_uniq30nt','UB40e1_filtered_uniq30nt','YA00e1_filtered_uniq30nt','YA05e1_filtered_uniq30nt','YA10e1_filtered_uniq30nt','YA20e1_filtered_uniq30nt','YA40e1_filtered_uniq30nt','YB00e1_filtered_uniq30nt','YB05e1_filtered_uniq30nt','YB10e1_filtered_uniq30nt','YB20e1_filtered_uniq30nt','YB40e1_filtered_uniq30nt'], axis=1)

    ''' decay2 '''

    file5 = open('/var/kristoph_flask/data/decay2_cumPa.txt', 'r')

    for line5 in file5:
        if re.search(gene_from_genome, line5):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneData.txt', 'a+') as f:
                f.write(line5)

    file6 = open('/var/kristoph_flask/data/decay2_paProb.txt', 'r')

    for line6 in file6:
        if re.search(gene_from_genome, line6):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneProb.txt', 'a+') as fa:
                fa.write(line6)

    data3 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneData.txt', sep = '\s+', header = None, names = ['#gene', 'chromo','position','strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e2_filtered_uniq30nt', 'UA05e2_filtered_uniq30nt', 'UA10e2_filtered_uniq30nt', 'UA20e2_filtered_uniq30nt', 'UA40e2_filtered_uniq30nt', 'UB00e2_filtered_uniq30nt', 'UB05e2_filtered_uniq30nt', 'UB10e2_filtered_uniq30nt', 'UB20e2_filtered_uniq30nt', 'UB40e2_filtered_uniq30nt', 'YA00e2_filtered_uniq30nt', 'YA05e2_filtered_uniq30nt','YA10e2_filtered_uniq30nt','YA20e2_filtered_uniq30nt','YA40e2_filtered_uniq30nt','YB00e2_filtered_uniq30nt','YB05e2_filtered_uniq30nt','YB10e2_filtered_uniq30nt','YB20e2_filtered_uniq30nt','YB40e2_filtered_uniq30nt','all'])

    df3 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneProb.txt', sep = '\s+', header = None, names = ['#gene', 'chromo','position','strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e2_filtered_uniq30nt', 'UA05e2_filtered_uniq30nt', 'UA10e2_filtered_uniq30nt', 'UA20e2_filtered_uniq30nt', 'UA40e2_filtered_uniq30nt', 'UB00e2_filtered_uniq30nt', 'UB05e2_filtered_uniq30nt', 'UB10e2_filtered_uniq30nt', 'UB20e2_filtered_uniq30nt', 'UB40e2_filtered_uniq30nt', 'YA00e2_filtered_uniq30nt', 'YA05e2_filtered_uniq30nt','YA10e2_filtered_uniq30nt','YA20e2_filtered_uniq30nt','YA40e2_filtered_uniq30nt','YB00e2_filtered_uniq30nt','YB05e2_filtered_uniq30nt','YB10e2_filtered_uniq30nt','YB20e2_filtered_uniq30nt','YB40e2_filtered_uniq30nt','all'])

    data3 = data3.drop(['strand','avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e2_filtered_uniq30nt', 'UA05e2_filtered_uniq30nt', 'UA10e2_filtered_uniq30nt', 'UA20e2_filtered_uniq30nt', 'UA40e2_filtered_uniq30nt', 'UB00e2_filtered_uniq30nt', 'UB05e2_filtered_uniq30nt', 'UB10e2_filtered_uniq30nt', 'UB20e2_filtered_uniq30nt', 'UB40e2_filtered_uniq30nt', 'YA00e2_filtered_uniq30nt', 'YA05e2_filtered_uniq30nt','YA10e2_filtered_uniq30nt','YA20e2_filtered_uniq30nt','YA40e2_filtered_uniq30nt','YB00e2_filtered_uniq30nt','YB05e2_filtered_uniq30nt','YB10e2_filtered_uniq30nt','YB20e2_filtered_uniq30nt','YB40e2_filtered_uniq30nt'], axis = 1)

    df3 = df3.drop(['strand', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'UA00e2_filtered_uniq30nt', 'UA05e2_filtered_uniq30nt', 'UA10e2_filtered_uniq30nt', 'UA20e2_filtered_uniq30nt', 'UA40e2_filtered_uniq30nt', 'UB00e2_filtered_uniq30nt', 'UB05e2_filtered_uniq30nt', 'UB10e2_filtered_uniq30nt', 'UB20e2_filtered_uniq30nt', 'UB40e2_filtered_uniq30nt', 'YA00e2_filtered_uniq30nt', 'YA05e2_filtered_uniq30nt','YA10e2_filtered_uniq30nt','YA20e2_filtered_uniq30nt','YA40e2_filtered_uniq30nt','YB00e2_filtered_uniq30nt','YB05e2_filtered_uniq30nt','YB10e2_filtered_uniq30nt','YB20e2_filtered_uniq30nt','YB40e2_filtered_uniq30nt'], axis = 1)

    ''' Steinmetz '''

    file7 = open('/var/kristoph_flask/data/steinmetz_cumPa.txt', 'r')

    for line7 in file7:
        if re.search(gene_from_genome, line7):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneData.txt', 'a+') as f:
                f.write(line7)

    file8 = open('/var/kristoph_flask/data/steinmetz_paProb.txt', 'r')

    for line8 in file8:
        if re.search(gene_from_genome, line8):
            with open('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneProb.txt', 'a+') as fa:
                fa.write(line8)

    data4 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneData.txt', sep = '\s+', header = None, names = ['#gene', 'chromo','position', 'strand', 'distToCDSstop','distToCDSstart', 'avgGcount', 'AGcount', 'Acount','Gcount', 'Lane8ByIB_30nt_Ttrim', 'Lane8ByIIB_30nt_Ttrim', 'Lane8ByIIIB_30nt_Ttrim', 'Lane8TS1248IB_30nt_Ttrim', 'Lane8TS1248IIIB_30nt_Ttrim', 'Lane8TS685IB_30nt_Ttrim', 'Lane8TS685IIB_30nt_Ttrim', 'Lane8TS801IB_30nt_Ttrim', 'Lane8TS801IIB_30nt_Ttrim', 'Lane8TS801IIIB_30nt_Ttrim', 'lane3BY4741III_30nt_Ttrim', 'lane3TS1248II_30nt_Ttrim', 'lane3TS685III_30nt_Ttrim', 'lane3TS685II_30nt_Ttrim', 'lane3TS801III_30nt_Ttrim', 'lane3TS801I_30nt_Ttrim', 'all'])

    df4 = pd.read_csv('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneProb.txt', sep = '\s+', header = None, names = ['#gene', 'chromo','position','strand', 'distToCDSstop', 'distToCDSstart', 'avgGcount', 'AGcount', 'Acount','Gcount', 'Lane8ByIB_30nt_Ttrim', 'Lane8ByIIB_30nt_Ttrim', 'Lane8ByIIIB_30nt_Ttrim', 'Lane8TS1248IB_30nt_Ttrim', 'Lane8TS1248IIIB_30nt_Ttrim', 'Lane8TS685IB_30nt_Ttrim', 'Lane8TS685IIB_30nt_Ttrim', 'Lane8TS801IB_30nt_Ttrim', 'Lane8TS801IIB_30nt_Ttrim', 'Lane8TS801IIIB_30nt_Ttrim', 'lane3BY4741III_30nt_Ttrim', 'lane3TS1248II_30nt_Ttrim', 'lane3TS685III_30nt_Ttrim', 'lane3TS685II_30nt_Ttrim', 'lane3TS801III_30nt_Ttrim', 'lane3TS801I_30nt_Ttrim', 'all'])

    data4  = data4.drop(['strand', 'avgGcount', 'AGcount', 'Acount', 'Gcount', 'Lane8ByIB_30nt_Ttrim', 'Lane8ByIIB_30nt_Ttrim', 'Lane8ByIIIB_30nt_Ttrim', 'Lane8TS1248IB_30nt_Ttrim', 'Lane8TS1248IIIB_30nt_Ttrim', 'Lane8TS685IB_30nt_Ttrim', 'Lane8TS685IIB_30nt_Ttrim', 'Lane8TS801IB_30nt_Ttrim', 'Lane8TS801IIB_30nt_Ttrim', 'Lane8TS801IIIB_30nt_Ttrim', 'lane3BY4741III_30nt_Ttrim', 'lane3TS1248II_30nt_Ttrim', 'lane3TS685III_30nt_Ttrim', 'lane3TS685II_30nt_Ttrim', 'lane3TS801III_30nt_Ttrim', 'lane3TS801I_30nt_Ttrim'], axis = 1)

    df4 = df4.drop(['strand', 'avgGcount', 'AGcount', 'Acount','Gcount', 'Lane8ByIB_30nt_Ttrim', 'Lane8ByIIB_30nt_Ttrim', 'Lane8ByIIIB_30nt_Ttrim', 'Lane8TS1248IB_30nt_Ttrim', 'Lane8TS1248IIIB_30nt_Ttrim', 'Lane8TS685IB_30nt_Ttrim', 'Lane8TS685IIB_30nt_Ttrim', 'Lane8TS801IB_30nt_Ttrim', 'Lane8TS801IIB_30nt_Ttrim', 'Lane8TS801IIIB_30nt_Ttrim', 'lane3BY4741III_30nt_Ttrim', 'lane3TS1248II_30nt_Ttrim', 'lane3TS685III_30nt_Ttrim', 'lane3TS685II_30nt_Ttrim', 'lane3TS801III_30nt_Ttrim', 'lane3TS801I_30nt_Ttrim'], axis = 1)


    d1 = []
    for m in range(len(data)):
        if (m == 0):
            d1.append((data['all'][m]))
        else:
            d1.append((data['all'][m]-data['all'][m-1]))

    data['p_indep'] = d1

    d2 = []
    for h in range(len(data2)):
        if (h == 0):
            d2.append((data2['all'][h]))
        else:
            d2.append((data2['all'][h]-data2['all'][h-1]))

    data2['p_indep'] = d2

    d3 = []
    for r in range(len(data3)):
        if (r == 0):
            d3.append((data3['all'][r]))
        else:
            d3.append((data3['all'][r]-data3['all'][r-1]))

    data3['p_indep'] = d3

    d4 = []
    for w in range(len(data4)):
        if (w == 0):
            d4.append((data4['all'][w]))
        else:
            d4.append((data4['all'][w]-data4['all'][w-1]))

    data4['p_indep'] = d4

    if (data.position[0] > data.position[1]):
        antisense = True
    else:
        antisense = False

    print(type(upstreamBuf))
    upstreamBuf = int(upstreamBuf)
    print(type(upstreamBuf))

    print(df_comp)

    if (antisense==True):
        print("yoyoyo")
        CDS1 = data['position'][0] + data['distToCDSstart'][0]
        print("nsync")
        gen_pos = CDS1 + upstreamBuf
        print("yost")
        print(type(df_comp.Position))
        df_comp['aligned_pos'] = gen_pos - df_comp.Position
        print("fak")
    elif (antisense==False):
        print("hihihi")
        CDS1 = data['position'][0] - data['distToCDSstart'][0]
        gen_pos = CDS1 - upstreamBuf
        df_comp['aligned_pos'] = df_comp.Position + gen_pos

    print(antisense)

    cumul2 = []

    for y in range(len(df_comp)):
        if (y == 0):
            print("quo")
            cumul2.append(df_comp.pASiteT[0])
        else:
            cumul2.append(cumul2[y-1]+(1 - cumul2[y-1])*df_comp.pASiteT[y])

    print("yeet")
    df_comp['cumul2'] = cumul2
    df_comp['cumul2'] = df_comp['cumul2']/df_comp.cumul2[len(df_comp)-1]

    if (antisense == True):
        print("bye")
        cds1 = data.position[0]+data.distToCDSstart[0]
        cds2 = data.position[0]+data.distToCDSstop[0]
    elif (antisense == False):
        print("poop")
        cds1 = data.position[0]-data.distToCDSstart[0]
        cds2 = data.position[0]-data.distToCDSstop[0]
    
    if (antisense == True):
        indep_c = figure(x_range = (df_comp['aligned_pos'][0], df_comp['aligned_pos'][len(df_comp)-1]))
        c0 = indep_c.line(data['position'], data['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        c1 = indep_c.line(data2['position'], data2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        c2 = indep_c.line(data3['position'], data3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        c3 = indep_c.line(data4['position'], data4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        c4 = indep_c.line(df_comp['aligned_pos'], df_comp['cumulative'], color = '#000000', legend = '?', alpha = 1)
        indep_c.add_layout(Arrow(end=VeeHead(size=25, fill_alpha=0.4, line_alpha=0), line_color="black", x_start=cds1, y_start=-0.05, x_end=cds2, y_end=-0.05, line_alpha = 0.4, line_width = 20))
        indep_c.below[0].formatter.use_scientific = False
        indep_c.legend.location = "top_left"
        indep_c.legend.click_policy="hide"

        indep_p = figure(x_range = indep_c.x_range)
        p0 = indep_p.line(data['position'], data['p_indep'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        p1 = indep_p.line(data2['position'], data2['p_indep'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        p2 = indep_p.line(data3['position'], data3['p_indep'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        p3 = indep_p.line(data4['position'], data4['p_indep'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        p4 = indep_p.line(df_comp['aligned_pos'], df_comp['pASiteT'], color = '#000000', legend = '?', alpha = 1)
        indep_p.below[0].formatter.use_scientific = False
        indep_p.legend.location = "top_left"
        indep_p.legend.click_policy="hide"

        staged_c = figure(x_range = (df_comp['aligned_pos'][0], df_comp['aligned_pos'][len(df_comp)-1]))
        sc0 = staged_c.line(data['position'], data['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        sc1 = staged_c.line(data2['position'], data2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        sc2 = staged_c.line(data3['position'], data3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        sc3 = staged_c.line(data4['position'], data4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        sc4 = staged_c.line(df_comp['aligned_pos'], df_comp['cumul2'], color = '#000000', legend = '?', alpha = 1)
        staged_c.add_layout(Arrow(end=VeeHead(size=25, fill_alpha=0.4, line_alpha=0), line_color="black", x_start=cds1, y_start=-0.05, x_end=cds2, y_end=-0.05, line_alpha = 0.4, line_width = 20))
        staged_c.below[0].formatter.use_scientific = False
        staged_c.legend.location = "top_left"
        staged_c.legend.click_policy="hide"

        staged_p = figure(x_range = staged_c.x_range)
        sp0 = staged_p.line(df['position'], df['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        sp1 = staged_p.line(df2['position'], df2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        sp2 = staged_p.line(df3['position'], df3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        sp3 = staged_p.line(df4['position'], df4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        sp4 = staged_p.line(df_comp['aligned_pos'], df_comp['pASiteT'], color = '#000000', legend = '?', alpha = 1)
        staged_p.below[0].formatter.use_scientific = False
        staged_p.legend.location = "top_left"
        staged_p.legend.click_policy="hide"

    elif (antisense == False):
        indep_c = figure()
        c0 = indep_c.line(data['position'], data['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        c1 = indep_c.line(data2['position'], data2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        c2 = indep_c.line(data3['position'], data3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        c3 = indep_c.line(data4['position'], data4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        c4 = indep_c.line(df_comp['aligned_pos'], df_comp['cumulative'], color = '#000000', legend = '?', alpha = 1)
        indep_c.add_layout(Arrow(end=VeeHead(size=25, fill_alpha=0.4, line_alpha = 0), line_color="black", x_start=cds1, y_start=-0.05, x_end=cds2, y_end=-0.05, line_alpha = 0.4, line_width = 20))
        indep_c.below[0].formatter.use_scientific = False
        indep_c.legend.location = "top_left"
        indep_c.legend.click_policy="hide"

        indep_p = figure(x_range = indep_c.x_range)
        p0 = indep_p.line(data['position'], data['p_indep'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        p1 = indep_p.line(data2['position'], data2['p_indep'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        p2 = indep_p.line(data3['position'], data3['p_indep'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        p3 = indep_p.line(data4['position'], data4['p_indep'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        p4 = indep_p.line(df_comp['aligned_pos'], df_comp['pASiteT'], color = '#000000', legend = '?', alpha = 1)
        indep_p.below[0].formatter.use_scientific = False
        indep_p.legend.location = "top_left"
        indep_p.legend.click_policy="hide"

        staged_c = figure()
        sc0 = staged_c.line(data['position'], data['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        sc1 = staged_c.line(data2['position'], data2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        sc2 = staged_c.line(data3['position'], data3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        sc3 = staged_c.line(data4['position'], data4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        sc4 = staged_c.line(df_comp['aligned_pos'], df_comp['cumul2'], color = '#000000', legend = '?', alpha = 1)
        staged_c.add_layout(Arrow(end=VeeHead(size=25, fill_alpha=0.4, line_alpha=0), line_color="black", x_start=cds1, y_start=-0.05, x_end=cds2, y_end=-0.05, line_alpha = 0.4, line_width = 20))
        staged_c.below[0].formatter.use_scientific = False
        staged_c.legend.location = "top_left"
        staged_c.legend.click_policy="hide"

        staged_p = figure()
        sp0 = staged_p.line(df['position'], df['all'], color = '#4dac26', legend = 'pcf11_DRS', alpha = 1)
        sp1 = staged_p.line(df2['position'], df2['all'], color = '#b8e186', legend = 'Decay1', alpha = 1)
        sp2 = staged_p.line(df3['position'], df3['all'], color = '#f1b6da', legend = 'Decay2', alpha = 1)
        sp3 = staged_p.line(df4['position'], df4['all'], color = '#d01c8b', legend = 'Steinmetz', alpha = 1)
        sp4 = staged_p.line(df_comp['aligned_pos'], df_comp['pASiteT'], color = '#000000', legend = '?', alpha = 1)
        staged_p.below[0].formatter.use_scientific = False
        staged_p.legend.location = "top_left"
        staged_p.legend.click_policy="hide"

    print("shoop")
    col1 = column(children=[indep_c, indep_p], sizing_mode='stretch_both')
    col2 = column(children=[staged_c, staged_p], sizing_mode='stretch_both')
    col3 = column(children=[s1, s2, s3, s4, s5], sizing_mode='stretch_both')

    print(os.listdir('/var/kristoph_flask/data/'))

    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneData.txt')
    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_pcf11_specGeneProb.txt')

    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneData.txt')
    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_decay_specGeneProb.txt')

    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneData.txt')
    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_decay2_specGeneProb.txt')

    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneData.txt')
    os.remove('/var/kristoph_flask/data/'+gene_from_genome+'_steinmetz_specGeneProb.txt')

    return(col1, col2, col3)


if __name__ == '__main__':
    application.run(debug=True)
