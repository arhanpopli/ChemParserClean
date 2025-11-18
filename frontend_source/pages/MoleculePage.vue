<template>
  <q-page padding>

    <div v-if="store.showSpinner">
      <span class="absolute-center">
        <q-spinner 
          color="primary"
          size="3em"
        />
      </span>
    </div>

    <div v-else>
      <MolErrorBanner />

      <div class="row q-mt-xl">

        <div
          class="col-12 col-lg"
          style="margin: auto 0;"
        >
          <Composer />

        </div>

        <div
          class="col-12 col-lg"
        >
          <UserInput />
        </div>

        <div
          class="col-12 col-lg"
         >
           <PDFViewer />
        </div>

      </div>
    </div>

  </q-page>
</template>

<script setup>
import { onBeforeMount, onMounted } from 'vue'
import Composer from 'components/Molecule/Composer.vue'
import UserInput from 'components/Molecule/UserInput.vue'
import PDFViewer from 'components/Molecule/PDFViewer.vue'
import MolErrorBanner from 'components/Molecule/ErrorBanner.vue'

import { useMoleculeStore } from 'src/stores/molecule'

const store = useMoleculeStore()


onBeforeMount (() => {
  const store = useMoleculeStore()
  store.showSpinner= false
  console.log(store.showSpinner)
})

// onMounted(() => {
//   store.showSpinner= false
// })

</script>
