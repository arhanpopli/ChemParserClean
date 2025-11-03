<template>
  <div>
    <q-dialog
      v-model="props.store.dialog"
      persistent
      :maximized="true"
      transition-show="slide-up"
      transition-hide="slide-down"
    >
      <q-card>
        <q-bar>
          <div
            class="col text-center">
            Press close button when you are done
          </div>

          <q-space />

          <q-btn
            icon="close"
            color="red"
            @click='props.store.closeComposer'
          />
        </q-bar>

        <q-card-section>
          <div
            :id="props.store.maxComposerId"
            :style="props.store.maxComposerStyle"
          >
          </div>

          <div class="text-center">
            <ConvertButton
              label="Convert"
              color="primary"
              @click=props.store.convert(props.store.maxComposerId)
            />
          </div>

        </q-card-section>

      </q-card>
    </q-dialog>
  </div>

</template>

<script setup>
// Maximized composer

import { onMounted } from 'vue'
import { setComposer } from 'src/hooks/setComposer'
import { useMoleculeStore } from 'src/stores/molecule'
import ConvertButton from 'src/components/Shared/Button.vue'

/*const molStore = useMoleculeStore()*/
const props = defineProps(['store'])


onMounted(() => {
  setTimeout(() => {
    setComposer(
      props.store.maxComposerId,
      props.store.maxComposerType,
      props.store.maxComposerData
    )
  }, 2)
})

</script>
