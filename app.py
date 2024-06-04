from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import dash_mantine_components as dmc
from analysis import event_table, submit_jobs, process_output
import zipfile
from flask import send_file

grnas = pd.read_csv('demo_CDK2.txt', delimiter="\t")

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(external_stylesheets=external_stylesheets)

marks = {}
for i in range(15):
    marks[i] = str(i)
e_df = event_table(grnas)
e_df = e_df[e_df["mutation_category"] != "Silent"]
fig_1 = px.strip(e_df,x="edit_location", y="pam_pattern",color="new_aa",stripmode = "overlay", template="plotly_white",
                     hover_data = ["sgRNA sequence","Edit","PAM"],
                     labels={
                         "edit_location":"Amino acid position",
                         "new_aa" : "New amino acid",
                         "pam_pattern": "PAM"
                     })
fig_2 = px.density_heatmap(e_df,
                               x="original_aa", 
                               y="new_aa", 
                               color_continuous_scale ="Blues",
                               template="plotly_white",
                               labels = {
                                   "original_aa":"Original amino acid",
                                   "new_aa":"Edited amino acid"
                               })
    
def format_list(ls):
    return [{"label": str(i), "value": i} for i in ls]

app.layout = dmc.MantineProvider(
    theme={
        "fontFamily": "'Inter', sans-serif",
        "primaryColor": "indigo",
        "components": {
            "Button": {"styles": {"root": {"fontWeight": 400}}},
            "Alert": {"styles": {"title": {"fontWeight": 500}}},
            "AvatarGroup": {"styles": {"truncated": {"fontWeight": 500}}},
        },
    },
    children =
    html.Div([dmc.Center(style={"height": 100, "width": "100%"},
                         children=dmc.Title("Base editor gRNA design", order=1),
                         ),
        dmc.Grid(children=[
                dmc.Col(span="auto", children = [
                    dmc.Stack(children=[
                        dmc.TextInput(label="Gene symbol",
                                        id="gene-symbol",
                                        type="text",
                                        value="CDK2"),
                        dmc.TextInput(label="ENSEMBL transcript id",
                                        id="transcript-name", 
                                        type="text", 
                                        value="ENST00000266970"),
                        dmc.CheckboxGroup(
                            id="pam-checklist",
                            label="Select protospace adjacent motif (PAM) types", 
                            orientation = "horizontal", 
                            children=[
                                dmc.Checkbox(label="NGG", value="NGG"),
                                dmc.Checkbox(label="NG", value="NG"),
                                dmc.Checkbox(label="NNN", value="NNN")
                            ],
                            value=["NGG", "NG"]),
                    ])
                ]),
                dmc.Col(span="auto", children = [
                    dmc.Stack(children=[ 
                        dmc.NumberInput(
                            id="intron-buffer",
                            label="Intron buffer",
                            value=20,
                            min=0,
                            step=1,
                        ),
                        dmc.NumberInput(
                            id="sg-len",
                            label="sgRNA length",
                            value=20,
                            min=1,
                            step=1,
                        ),
                        dmc.CheckboxGroup(
                            label= "Filter edits in GC motif regions",
                            id="filter-gc-content",
                            orientation= "horizontal",
                            children=dmc.Checkbox(
                                label="Filter GC content",
                                value="Filter GC content",
                                mb=10),
                            value = []
                        )
                    ])
                ]),
            dmc.Col(span="auto",children =[
                dmc.Stack([
                    dmc.CheckboxGroup(
                        label="Select edit types",
                        orientation="horizontal",
                        id="edit-types",
                        children=[
                            dmc.Checkbox(
                                label="A-T",
                                value="A-T",
                                mb=10
                            ),
                            dmc.Checkbox(
                                label="G-C",
                                value="G-C",
                                mb=10
                            )
                        ],
                        value=["A-T","G-C"]
                    ),
                    dmc.Text("Edit window (nucleotides)"),
                    dcc.RangeSlider(
                        id="range-slider-callback",
                        pushable=2,
                        step=1,
                        value=[4,8],
                        min=0,
                        max=15,
                        marks = marks,
                    ),
                dmc.Button("Submit",
                           id="submit-button", 
                           variant="gradient", 
                           loaderProps={"type": "dots"}),
                ]),
            ]),
        ]),
        dmc.Space(h=10),
        dmc.Stack([
            dmc.Divider(variant="dashed"),
            dmc.Grid(children=[
                dmc.Col(span="auto", children = dcc.Graph(id='dot-plot',figure=fig_1)),
                dmc.Col(span=3, children = dcc.Graph(id="heat-map",figure=fig_2))
            ]),
            dmc.Button("Download",
                       id="download-button",
                       variant="gradient")
        ])
    ])
)

@app.callback(
    Output('dot-plot', 'figure'),
    Output('heat-map', 'figure'),
    Input('submit-button','n_clicks'),
    State("pam-checklist",'value'),
    State("edit-types",'value'),
    State("intron-buffer",'value'),
    State("filter-gc-content",'value'),
    State("range-slider-callback",'value'),
    State("transcript-name",'value'),
    State("gene-symbol",'value'),
    State("sg-len",'value'),
    prevent_initial_call=True,
    running=[(Output("submit-button", "loading"), True, False)]
)
def submit_button_click(clicks,
                        pam_types, 
                        edit_types,
                        intron_buffer, 
                        filter_gc_content, 
                        edit_window,
                        transcript_id,
                        gene_symbol,
                        sg_len):
    if clicks is not None:
        fil_gc = len(filter_gc_content) != 0
        d = submit_jobs(pam_types,edit_types,intron_buffer,fil_gc,"-".join(str(x) for x in edit_window),transcript_id,gene_symbol,sg_len)
        df = process_output(d)
        with zipfile.ZipFile('grna_results.zip', 'w') as myzip:
            for f in d.keys():   
                myzip.write(f)
        e_df = event_table(df)
        e_df = e_df[e_df["mutation_category"] != "Silent"]
        fig_1 = px.strip(e_df,x="edit_location", y="pam_pattern",color="new_aa",stripmode = "overlay", template="plotly_white",
                            hover_data = ["sgRNA sequence","Edit","PAM"],
                            labels={
                                "edit_location":"Amino acid position",
                                "new_aa" : "New amino acid",
                                "pam_pattern": "PAM"
                            })
        fig_2 = px.density_heatmap(e_df,
                                    x="original_aa", 
                                    y="new_aa", 
                                    color_continuous_scale ="Blues",
                                    template="plotly_white",
                                    labels = {
                                        "original_aa":"Original amino acid",
                                        "new_aa":"Edited amino acid"
                                    })
        return(fig_1,fig_2)

@app.callback(
    Output("download-button","value"),
    Input("download-button","n_clicks"),
    prevent_initial_call=True
)
def download():
    send_file( filename=
        "./grna_results.zip")
    return "Download"

if __name__ == '__main__':
    app.run_server(debug=True)