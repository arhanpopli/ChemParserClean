<template>
  <div
    :id="props.store.composerId"
    :style="props.store.composerStyle"
  >

  </div>
  <div class="text-center q-mt-md">
    <ConvertButton
      label="Convert"
      color="primary"
      @click=props.store.convert(props.store.composerId)
      :loading="props.store.isConverting"
    />
  </div>
</template>

<script setup>

import { onMounted, onBeforeUnmount } from 'vue'
import { useMoleculeStore } from 'src/stores/molecule'
import { setComposer } from 'src/hooks/setComposer'
import ConvertButton from 'src/components/Shared/Button.vue'

const props = defineProps(["store"])

onMounted(() => {

  setComposer(props.store.composerId,
              props.store.composerType,
              props.store.composerData
             )
})

// Disabled while working on reaction. Needs to be back.
onBeforeUnmount(() => {
      props.store.beforeUnmount()
    })

</script>
