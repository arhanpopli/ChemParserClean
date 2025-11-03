import { defineStore } from 'pinia';
import isNil from 'lodash/isNil'
import isEmpty from 'lodash/isEmpty'
import { getDocObject, loadDocToComposer  } from 'src/hooks/getDocObject.js'
import { copyData } from 'src/hooks/copyToClipboard.js'
import { useMakeRequest } from 'src/hooks/makeRequest.js'
import { isRpresent } from 'src/hooks/checkRsubgroup.js'
import { useSettings } from 'src/stores/settings'

const settings = useSettings()

export const useReactionStore = defineStore('reaction', {
  state: () => ({
    composerId: "composer",
    composerType: "reaction",
    composerStyle: "width:100%;height:45vh;z-index: 1",
    composerData: null,

    // Maximized composer
    dialog: false,
    maxComposerId: "maxComposer",
    maxComposerType: "reaction",
    maxComposerStyle: "width:100%;height:90vh;z-index: 1",
    maxComposerData: null,

    // Text Area
    textAreaData: "",

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

    // Tooltip
    showTooltip: false,

    error: null,

    isConverting: false,

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

    getDocJsonAndMolFiles(compId) {

      let composer;

      if (compId === "composer"){
        composer = Kekule.Widget.getWidgetById(this.composerId)
      }

      else {
        composer = Kekule.Widget.getWidgetById(this.maxComposerId)
        this.closeComposer()
      }

      const docJSON = getDocObject(compId)
      const mols = composer.exportObjs(Kekule.Molecule)
      let molFiles = {}
      let is_R_present = false

      // Mapping molecule id (from kekule-json format) to MOL Format of the mol
      for (var i = 0; i < mols.length; ++i) {
        var mol = mols[i];

        let jsonMol = JSON.parse(Kekule.IO.saveFormatData(mol, "Kekule-JSON"))
        // Formula tool creates "Molecule type" as well (C6H5OH) but MOL File
        // cannot be generated.
        if (jsonMol.hasOwnProperty("formula")) {
          continue
        }
        // Check if R (substituent group) is present since mol2chemfig
        // cannot parse it
        let MOLFile = Kekule.IO.saveFormatData(mol, 'mol')
        // Check if R (substituent group) is present in the structure
        if (isRpresent(MOLFile)) {
          is_R_present = true
        }

        molFiles[jsonMol['id']] = MOLFile
      }

      let payload = {
        "docJSON": docJSON,
        "mol_files": molFiles,
        "Rpresent": is_R_present
      }

      return payload
    },

    async convert(compId) {

      const payload = this.getDocJsonAndMolFiles(compId)

      if (payload.Rpresent) {
        this.error = "Mol2chemfig cannot convert R group. Rreplace it with other group"
        return null
      }

      const url_ = "/m2cf/reaction/convert"

      this.isConverting = true

      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if (isNil(response_data)) {
        this.isConverting = false
        this.error = "Something went wrong while converting to chemfig"
      }
      if (isNil(error)) {
        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
        this.error = response_data.error
        this.isConverting = false
      } else {
        this.error = error
        this.isConverting = false
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

      let payload = this.getDocJsonAndMolFiles('composer')

      const options = {
        "selections": this.selection,
        "angle": this.angle,
        "indentation": this.indentation,
        "h2": this.h2
      }

      payload = {...payload, ...options}

      const url_ = "/m2cf/reaction/apply"

      this.isApplying = true

      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if ( isNil(error) ) {
        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
        this.isApplying = false
        this.error = response_data.error

      } else {
        console.log(error)
        this.isApplying = false
        this.error = error
      }
    },

    async reset() {
      const payload = this.getDocJsonAndMolFiles('composer')
      const url_ = "/m2cf/reaction/reset"

      this.isResetting = true

      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if ( isNil(error) ) {

        this.textAreaData = response_data.chemfig
        this.src = response_data.pdflink
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

    async updateChemfig() {
      if (isEmpty(this.textAreaData)) {
        this.error = "Cannot submit empty data"
        return
      }
      const url_ = "/m2cf/reaction/update_chemfig"

      const payload = {
        "textAreaData": this.textAreaData
      }

      const { response_data, error, isLoading } = await useMakeRequest(url_, payload)
      if ( isNil(error) && isNil(response_data.error) ) {
        console.log(response_data)
        this.src = response_data.pdflink
      } else {
        this.error = response_data.error
      }
    },

    beforeUnmount(){
      const keepOnLeave = settings.filterByTerm('tag',
                                                   'keepOnLeave',
                                                   'reaction'
      )
      if ( keepOnLeave ) {
        this.composerData = getDocObject(this.composerId)
      }
      else {
        this.composerData = ""
      }
    }

  },
});
