<template>
  <div>
    <MaximizeButton :store="reactionStore"/>
    <Composer :store="reactionStore" />
    <MaxComposer
      v-if="reactionStore.dialog"
      :store="reactionStore"
    />
  </div>
</template>

<script setup>

import { useReactionStore } from 'src/stores/reaction'

import Composer from  'src/components/Shared/Composer.vue'
import MaxComposer from  'src/components/Shared/MaxComposer.vue'
import MaximizeButton from 'src/components/Shared/MaximizeButton.vue'

const reactionStore = useReactionStore()

</script>
