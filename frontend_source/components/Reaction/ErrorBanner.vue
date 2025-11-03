<template>
  <ErrorBanner :store="reactionStore" />
</template>

<script setup>

import ErrorBanner from 'src/components/Shared/ErrorBanner.vue'
import { useReactionStore } from 'src/stores/reaction'

const reactionStore = useReactionStore()

</script>
