import { defineStore } from 'pinia';
import { LocalStorage } from 'quasar'

export const useSettings = defineStore('Settings', {
  state: () => ({
    settings: [

      {
        id_: 0,
        description: "Insert a searched chemical structure into Composer",
        type: "molecule",
        tag: "intoComposer",
        model: false
      },

      {
        id_: 1,
        description: "Preserve a drawn structure in the composer/ sketcher when leaving the page",
        type: "molecule",
        tag: "keepOnLeave",
        model: false
      },

      {
        id_: 2,
        description: "Hide tooltips for checkbox options",
        type: "molecule",
        tag: "hideTooltips",
        model: false
      },

      {
        id_: 3,
        description: "Preserve a drawn structure in the composer/ sketcher when leaving the page",
        type: "reaction",
        tag: "keepOnLeave",
        model: false
      },

      {
        id_: 4,
        description: "Hide tooltips for checkbox options",
        type: "reaction",
        tag: "hideTooltips",
        model: false
      }
    ],

  }),

  getters: {

  },

  actions: {

    setValue(item){
      this.settings[item.id_].model = item.model
      LocalStorage.set(item.id_, item.model)
    },

    filterByTerm(key, value, type){
      let filteredSettings = this.settings.filter(
        item => item[key] == value && item['type'] == type )
      return filteredSettings[0].model
    },

    getSettings(key, value){
      const settings = methods.filterByTerm(key, value, type)
      return settings
    },

    setSettings(){
      const localKeys = Object.keys(localStorage)
      const regex = /^[0-9]+$/
      if (localKeys.length > 0) {
        // try except continue if kekule saves data into localstorage
        try {
          for (let key of localKeys) {
            if (key.match(regex)) {
              const intKey = parseInt(key)
              let currentSetting = this.settings[intKey]
              const val = LocalStorage.getItem(key)
              currentSetting.model = val
            }
            else {
              continue
            }
          }
        }
        catch (error) {
          console.log(error)
        }
      }
    }
  }
})
