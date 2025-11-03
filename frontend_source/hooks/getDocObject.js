export const getDocObject = (composerId) => {
  const composer = Kekule.Widget.getWidgetById(composerId);
  const chemDoc = composer.getChemObj();
  return Kekule.IO.saveFormatData(chemDoc , 'Kekule-JSON');
}


export const loadDocToComposer = (composerId, composerData) => {
  const composer = Kekule.Widget.getWidgetById(composerId);
  const chemDoc = composer.getChemObj();
  const doc = Kekule.IO.loadFormatData(composerData, 'Kekule-JSON');
  composer.setChemObj(doc);
}

export const loadMolBlockToComposer = (composerId, molBlock) => {
  const composer = Kekule.Widget.getWidgetById(composerId);
  const doc = Kekule.IO.loadFormatData(molBlock, 'mol');
  composer.setChemObj(doc);

}
