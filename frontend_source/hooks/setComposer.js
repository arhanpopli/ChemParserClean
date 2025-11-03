const setComposer = (composerId, composerType, composerData) => {

  const composer = new Kekule.Editor.Composer(document.getElementById(composerId));
  composer.setCommonToolButtons(
    [
      'newDoc', 'loadData', 'saveData',
      'undo', 'redo', 'copy', 'cut',
      'paste', 'zoomIn', 'reset', 'zoomOut'
    ]
  );
  if ( composerType == 'molecule' ) {
    composer.setChemToolButtons(
      [
        'manipulate', 'erase', 'bond',
        'atom', 'ring', 'charge'
      ]
    )
  }
  composer.getRenderConfigs().getLengthConfigs().setDefBondLength(19);
  composer.getRenderConfigs().getLengthConfigs().setAtomFontSize(10);
  composer.getRenderConfigs().getColorConfigs().setGlyphStrokeColor('#000000');  // set to black
  composer.getRenderConfigs().getColorConfigs().setGlyphFillColor('#000000');
  composer.newDoc();
  if (composerData) {
    const doc = Kekule.IO.loadFormatData(composerData, 'Kekule-JSON');
    composer.setChemObj(doc);
  }
}

export  {setComposer}
