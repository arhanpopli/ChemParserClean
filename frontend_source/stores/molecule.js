import { defineStore } from 'pinia';
import { getDocObject,
         loadDocToComposer,
         loadMolBlockToComposer
       } from 'src/hooks/getDocObject.js'
import { useSettings } from 'src/stores/settings'
import isNil from 'lodash/isNil'
import isEmpty from 'lodash/isEmpty'
import {useMakeRequest} from 'src/hooks/makeRequest'
import { copyData } from 'src/hooks/copyToClipboard.js'
import { isRpresent } from 'src/hooks/checkRsubgroup.js'

const settings = useSettings()

export const useMoleculeStore = defineStore('molecule', {
  state: () => ({
    showSpinner: true,
    composerId: "composer",
    composerType: "molecule",
    composerStyle: "width:100%;height:450px;z-index: 1",
    composerData: null,

    // Maximized composer
    dialog: false,
    maxComposerId: "maxComposer",
    maxComposerType: "molecule",
    maxComposerStyle: "width:100%;height:90vh;z-index: 1",
    maxComposerData: null,

    // User Input
    searchTerm: "",
    textAreaData: "",
    textAreaLastData: "",
    lastChemData: "",
    lastChemFormat: "",

    // Copy Tooltip
    showCopyTooltip: false,

    // Checkboxes tooltips
    hideOptionsTooltips: false,

    // Iframe
    src: "",

    options: [
      {
        label: 'compact view',
        value: '-z',
        tooltip: "A more compact chemfig format"
      },
      {
        label: 'fancy bonds',
        value: '-f',
        tooltip: 'Show nicer double and triple bonds'
      },
      {
        label: 'aromatic',
        value: '-o',
        tooltip: 'Show circle in aromatic compounds instead of double bonds'
      },
      {
        label: 'show carbon',
        value: '-c',
        tooltip: 'Show carbon atoms as elements'
      },
      {
        label: 'show methyl',
        value: '-m',
        tooltip: 'Show methyl group as elements'
      },
      {
        label: 'flip',
        value: '-p',
        tooltip: 'Flip the structure horizontally'
      },
      {
        label: 'flop',
        value: '-q',
        tooltip: 'Flip the structure vertically'
      },
      {
        label: 'atom-numbers',
        value: '-n',
        tooltip: "Asign number to atoms except hydrogens"
      }
    ],

    selection: [],

    angle: 0,
    indentation: 4,

    h2: "keep",
    h2_options: ["keep", "add", "delete"],

    error: null,

    isConverting: false,

    isSubmitting: false,

    isSearching: false,

    isApplying: false,

    isResetting: false

  }),

  getters: {
  },

  actions: {

    closeBanner() {
      this.error = null
    },

    closeComposer() {
      const composerData = getDocObject(this.maxComposerId)
      loadDocToComposer(this.composerId, composerData)
      this.dialog = false
    },

    maximizeComposer(){
      this.maxComposerData = getDocObject(this.composerId)
      this.dialog = true
    },

    hideTooltips(model) {
      if (model) {
        this.hideOptionsTooltips = true
      }
      else {
        this.hideOptionsTooltips = false
      }
    },

    async convert(compId) {

      let composer;

      if (compId === "composer"){
        composer = Kekule.Widget.getWidgetById(this.composerId)
      }

      else {
        composer = Kekule.Widget.getWidgetById(this.maxComposerId)
        this.closeComposer()
      }

      const mols = composer.exportObjs(Kekule.Molecule)
      const mol = mols[0]
      // If convert btn is clicked when composer is empty/no molecule
      if ( isNil(mol) ) {
        this.error = "Please draw a structure before convering it to chemfig"
        return null
      }

      const molBlock = Kekule.IO.saveFormatData(mol, 'mol')

      if (isRpresent(molBlock)) {
        this.error = "Mol2chemfig cannot convert R group. Rreplace it with another group"
        return
      }

      let payload = {
        "data": molBlock,
        "type": "mol_block"
      }

      // const url_ = "/m2cf/convert"
      const url_ = "/m2cf/convert"
      this.isConverting = true
      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      // Check if drawn structure in sketcher is valid
      if (!isEmpty(response_data.error)){
        this.error = response_data.error
        this.isConverting = false
        return null
      }

      if (isNil(error)) {
        this.textAreaData = response_data.chemfig
        this.textAreaLastData = this.textAreaData
        this.isConverting = response_data.isLoading
        this.src = response_data.pdflink
      }

      else {
        this.isConverting = response_data.isLoading
        this.error = response_data.error
      }
      this.lastChemData = molBlock
      this.lastChemFormat = "mol"
    },

    setDataType(data){

      let dataType = null
      let payload = {}

      if (data.includes('chemfig')) {
        dataType = 'chemfig'
        payload["data"] = data
        this.lastChemFormat = dataType
      }

      else if ( data.includes('END')) {
        dataType = 'mol'
        payload["data"] = data
        this.lastChemData = data
        this.lastChemFormat = dataType
      }

      else {
        dataType = "smiles"
        payload["data"] = data
        this.lastChemData = data
        this.lastChemFormat = dataType
      }

      payload["dataType"] = dataType
      return [dataType, payload]
    },

    async updateChemfig(payload) {
      const url_ = "/m2cf/update_chemfig"
      this.isSubmitting = true
      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if ( isNil(error) ) {
        this.isLoading = false
        this.textAreaLastData = this.textAreaData
        this.lastChemFormat = 'chemfig'
        this.src = response_data.pdflink
        this.isSubmitting = false
      } else {
        console.log(error)
        this.isSubmitting = false
        this.error = response_data.error
      }
    },

    async submit() {

      if (isEmpty(this.textAreaData)) {
        this.error = "You cannot submit empty text area!"
        return null
      }

      const payload = {
        "textAreaData": this.textAreaData
      }

      if (this.textAreaData.includes('chemfig')) {
        this.updateChemfig(payload)
        return null
      }


      this.isSubmitting = true
      const url_ = "/m2cf/submit"
      const { response_data, error, isLoading} = await useMakeRequest(url_,
                                                        payload)
      if ( isNil(error) ) {
        this.textAreaData = response_data.chemfig
        this.lastChemFormat = response_data.chem_format
        this.lastChemData = response_data.chem_data
        this.src = response_data.pdflink
        this.isSubmitting = isLoading
      } else {
        this.error = response_data.error
        this.isSubmitting = false
      }
    },

    async search() {
      if (isEmpty(this.searchTerm)) {
        alert("Please enter a compound's name")
      }

      let payload = {
        "searchTerm": this.searchTerm,
      }

      const url_ = "/m2cf/search"
      this.isSearching = true
      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if ( isNil(error) ) {
        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
        //this.isSubmitting = isLoading
        this.textAreaLastData = this.textAreaData
        this.lastChemData = response_data.smiles
        this.lastChemFormat = 'smiles'
        this.error = response_data.error
        this.searchTerm = ""
        this.isSearching = false

        const intoComposer = settings.filterByTerm('tag', "intoComposer", "molecule")
        if (intoComposer) {
          loadMolBlockToComposer(this.composerId, response_data.molblock)
        }

      } else {
        this.isSearching = false
        console.log(error)
        this.error = error
      }
    },

    copyToClipBoard() {
      this.error = copyData(this.textAreaData)
      if (isNil(this.error)) {
        this.showCopyTooltip = true
        setTimeout(() => {
          this.showCopyTooltip = false
        }, 1000)
      }
    },

    async apply() {

      const url_ = "/m2cf/apply"

      const payload = {
        "chem_data": this.lastChemData,
        "chem_format": this.lastChemFormat,
        "selections": this.selection,
        "angle": this.angle,
        "indentation": this.indentation,
        "h2": this.h2
      }

      this.isApplying = true

      const { response_data, error, isLoading } = await useMakeRequest(url_,
      payload)

      if ( isNil(error) ) {
        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
        this.isApplying = false
        this.textAreaLastData = this.textAreaData
        this.error = response_data.error

      } else {
        console.log(error)
        this.isApplying = false
        this.error = error
      }
    },

    async reset() {


      const url_ = "/m2cf/reset"

      const payload = {
        "chem_data": this.lastChemData,
        "chem_format": this.lastChemFormat,
      }

      this.isResetting = true

      const { response_data, error, isLoading } = await useMakeRequest(url_,
      payload)

      if ( isNil(error) ) {
        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
        this.textAreaLastData = this.textAreaData
        this.error = response_data.error

        this.selection = []
        this.angle = 0
        this.indentation = 4
        this.h2 = "keep"

        this.isResetting = false

      } else {
        console.log(error)
        this.isResetting = false
        this.error = error
      }
    },

    beforeUnmount(){
      const keepOnLeave = settings.filterByTerm('tag', 'keepOnLeave', "molecule")
      if ( keepOnLeave ) {
        this.composerData = getDocObject(this.composerId)
      }
      else {
        this.composerData = ""
      }
    }
  },
});
