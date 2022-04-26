
//run: node Gene_T_JavaScipt.js
//Pre-Processing data
var np = require('numjs'); //numpy
var DataFrame=require('pandas-js');//pandas
var os = require('os'); //os 
var Plotly = require('plotly.js'); //Plotly for matplotlib

//Read TCGA dataset
let TCGA_Data = pd.read_csv('TCGA1.csv');
let idx_sample=TCGA_Data.sampleID;
TCGA_Data=TCGA_Data.drop('sampleID', axis=1);
let idx_Genes=TCGA_Data.columns;
TCGA_Data=pd.DataFrame(TCGA_Data.T.values,columns=idx_sample);
TCGA_Data.index=idx_Genes;

//Read Phenotype dataset
let TCGA_clinical=pd.read_csv('TCGA_clinical_metadata.csv');

//Read Layout 
let organic_layout= pd.read_csv('Layout.txt',header=None,sep='	');
organic_layout.index= organic_layout.iloc[":,0"];
organic_layout=organic_layout.drop(0, axis=1);
organic_layout.columns=['X','Y','PPI'];

//Merge TCGA Data and Layout
let Merge_dataset=pd.merge(organic_layout, TCGA_Data, left_index=True, right_index=True);

let standard_embedding= Merge_dataset[['X','Y']];
TCGA_Data=Merge_dataset[idx_sample];

//Create GeneTerrain for each patient
let res=200;
let x_ = np.linspace(0, 1, res);
let y_ = np.linspace(0, 1, res);
let X, Y = np.meshgrid(x_, y_);
let pi=3.14;

for (name_sample in TCGA_Data.columns){
      let gaussian=np.zeros((res,res));

      
      for (i in range(np.shape(standard_embedding)[0])){
          let sigma= .05;
          let x=standard_embedding['X'][i];
          let y=standard_embedding['Y'][i];
          let amp1=TCGA_Data[name_sample][i];
          
          let gaussian1 =np.exp(-(((X-x)/sigma)**2+((Y-y)/sigma)**2)); 
          let gaussian =amp1* gaussian1 + gaussian;}

      Plotly.Figure(figsize=(7, 6)); //Plotly.figure(figsize=(7, 6));
      norm=Plotly.Normalize(-10,10);
      Plotly.pcolormesh(agaussian, cmap = 'jet',norm=norm);
      Plotly.colorbar();
      idxxs=np.where(TCGA_clinical['sampleID']==name_sample)[0];
      category=TCGA_clinical['GeneExp_Subtype'][idxxs].tolist()[0];
 
      path=os.getcwd();
      
      if (category =='Classical'){
            Plotly.write_image(path+ '\\desktop\\GeneTerain_GBM\\Classical\\'+ name_sample +'.jpeg',dpi=100);
      }         
      else if (category=='Mesenchymal'){
            Plotly.write_image(path+'\\desktop\\GeneTerain_GBM\\Mesenchymal\\'+ name_sample +'.jpeg',dpi=100);
      }
      else if (category=='Neural'){
            Plotly.write_image(path+'\\desktop\\GeneTerain_GBM\\Neural\\'+ name_sample +'.jpeg', dpi=100);
      }
      else {
            Plotly.write_image(path+'\\desktop\\GeneTerain_GBM\\Proneural\\'+ name_sample +'.jpeg',dpi=100);
      }
      Plotly.show();     }
    
   
    
   