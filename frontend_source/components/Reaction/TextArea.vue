<template>
  <div>

    <TextArea
      class="full-width"
      rows='25'
      v-model="reactionStore.textAreaData"
      @keyup.ctrl.enter="reactionStore.updateChemfig"
    />

    <div class="float-right q-mr-lg">
      <q-btn
        flat
        round
        color="grey-8"
        icon="content_copy"
        @click="reactionStore.copyToClipBoard"
      >
      </q-btn>
    </div>

    <Tooltip :store="reactionStore"/>
  </div>
</template>

<script setup>

import TextArea from 'components/Shared/TextArea.vue'
import { useReactionStore } from 'src/stores/reaction'
import Tooltip from 'components/Shared/Tooltip.vue'

const reactionStore = useReactionStore()
</script>
