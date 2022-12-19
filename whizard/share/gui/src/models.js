export function SindarinModel(name, description) {
  this.modelName = name;
  this.description = description;
  this.toString = () => 'model = ' + this.modelName;
}


// TODO: (bcn 2016-06-30) parameters are not setable yet
export function SindarinModelData(modelString) {
  this.model = new SindarinModel(modelString);
  this.parameters = [];
  this.writeToSindarin = () => {
    let src = this.model.toString() + '\n';
    for (let j = 0; j < this.parameters.length; j++) {
      src += this.parameters.toString() + '\n';
    }
    return src;
  };
}


function fillModelList() {
  const modelList = [];
  const models = [
    {name: 'AltH', description: 'An SM extension for VV scattering'},
    {name: 'GravTest', description: 'SQED with gravitino'},
    {name: 'HSExt', description: '?????'},
    {name: 'Littlest', description: 'Littles Higgs'},
    {name: 'Littlest_Eta', description: 'Littlest Higgs with ungauged U(1)'},
    {name: 'Littlest_Tpar', description: 'Littlest Higgs with T parity'},
    {name: 'MSSM', description: 'Minimal supersymmetric standard model'},
    {name: 'MSSM_CKM', description: 'MSSM with CKM matrix'},
    {name: 'MSSM_Grav', description: 'MSSM with gravitinos'},
    {name: 'MSSM_Hgg', description: 'MSSM with Hgg-coupling (???)'},
    {name: 'NMSSM', description: 'Next-to-Minimal supersymmetric standard model'},
    {name: 'NMSSM_CKM', description: 'NMSSM with CKM matrix'},
    {name: 'NMSSM_Hgg', description: 'NMSSM with Hgg-coupling (???)'},
    {name: 'NoH_rx', description: 'SM with anomalous Higgs couplings'},
    {name: 'PSSSM', description: 'Extended SUSY models'},
    {name: 'QCD', description: 'QCD with d,u,s,c,b,t,g'},
    {name: 'QED', description: 'QED with e,μ,τ,γ'},
    {name: 'Simplest', description: 'Simplest Little Higgs (anomaly-free)'},
    {name: 'Simplest_univ', description: 'Simplest Little Higgs (universal)'},
    {name: 'SM', description: 'Standard model'},
    {name: 'SM_ac', description: 'SM with anomalous gauge couplings'},
    {name: 'SM_ac_CKM', description: 'SM with anomalous gauge couplings and CKM matrix'},
    {name: 'SM_CKM', description: 'SM with CKM matrix'},
    {name: 'SM_hadrons', description: '???'},
    {name: 'SM_Higgs', description: 'SM with Hgg, Hγγ and Hμμ'},
    {name: 'SM_rx', description: 'SM with anomalous Higgs couplings'},
    {name: 'SM_top', description: 'SM with chareg 4/3 top'},
    {name: 'SM_top_anom', description: 'SM with anomalous top coupings'},
    {name: 'SM_tt_threshold', description: 'SM with top-threshold resummation'},
    {name: 'SM_ul', description: 'SM with anomalous Higgs couplings'},
    {name: 'SSC', description: 'SM extension for VV scattering'},
    {name: 'SSC_2', description: 'SM extension for VV scattering'},
    {name: 'SSC_AltT', description: 'SM extension for VV scattering'},
    {name: 'Template', description: 'Augmentable SM template'},
    {name: 'THDM', description: 'Two-Higgs Doublet Model'},
    {name: 'THDM_CKM', description: 'Two-Higgs Doublet Model with CKM matrix'},
    {name: 'Threeshl', description: '???'},
    {name: 'Threeshl_nohf', description: '???'},
    {name: 'UED', description: 'Universal Extra Dimensions'},
    {name: 'Xdim', description: 'SM with graviton'},
    {name: 'Zprime', description: "SM with Z'"},
  ];
  for (const model of models) {
    modelList.push(new SindarinModel(model.name, model.description));
  }
  return modelList;
}

export function setupJquery(ToolbarColumns) {
  const Models = fillModelList();
  for (let k = 0; k < Models.length; k += ToolbarColumns) {
    $('#pop_models').append('<div class="row">');
    for (let i = k; i < k + ToolbarColumns; i++) {
      const modelName = (Models[i] === undefined) ? '&nbsp;' : Models[i].modelName;
      const modelDescription = (Models[i] === undefined) ?
        '&nbsp;' : Models[i].description;
      $('#pop_models').append('<div class="col-md-3"><a href="javascript:;" title="'
          + modelDescription + '" class="model">' + modelName + '</a></div>');
    }
    $('#pop_models').append('</div>');
  }
}
