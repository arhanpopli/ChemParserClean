<template>
  <Iframe
    :src="reactionStore.src"
    style="width:100%;height:33vh"
  />
</template>

<script setup>

import Iframe from 'components/Shared/Iframe.vue'
import { useReactionStore } from "src/stores/reaction"

const reactionStore = useReactionStore()

</script>
