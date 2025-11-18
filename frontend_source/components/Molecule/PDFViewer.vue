<template>
  <div class='q-mt-xl'>
    <div id="pdf_viewer" class="row justify-center">
      <Iframe
        :src="molStore.src"
        style="width:80%;height:35vh"
      />
    </div>

    <div>
      <Selections :store="molStore"/>
    </div>
   <div>
     <ApplyReset
       :apply="molStore.apply"
       :reset="molStore.reset"
       :isApplying = "molStore.isApplying"
       :isResetting = "molStore.isResetting"
     />
   </div>
  </div>
</template>

<script setup>

import { useMoleculeStore } from "src/stores/molecule"
import Iframe from 'components/Shared/Iframe.vue'
import Selections from 'components/Shared/Selections.vue'
import ApplyReset from 'components/Shared/ApplyReset.vue'

import { LocalStorage } from 'quasar'
import { useSettings } from 'src/stores/settings'

const molStore = useMoleculeStore()
const settings = useSettings()
const shouldHidetooltip = settings.filterByTerm("tag", "hideTooltips", "molecule")
molStore.hideTooltips(shouldHidetooltip)

</script>
