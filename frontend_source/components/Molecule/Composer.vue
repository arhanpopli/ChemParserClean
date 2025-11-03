<template>
  <div>
    <MaximizeButton :store="molStore"/>
    <Composer :store="molStore"/>
    <MaxComposer
      v-if="molStore.dialog"
      :store="molStore"
    />
  </div>
</template>

<script setup>

import { useMoleculeStore } from 'src/stores/molecule'

import Composer from  'src/components/Shared/Composer.vue'
import MaxComposer from  'src/components/Shared/MaxComposer.vue'
import MaximizeButton from 'src/components/Shared/MaximizeButton.vue'

const molStore = useMoleculeStore()

</script>
