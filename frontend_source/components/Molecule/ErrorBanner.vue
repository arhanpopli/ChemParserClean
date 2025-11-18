<template>
  <ErrorBanner :store="molStore" />
</template>

<script setup>

import ErrorBanner from 'src/components/Shared/ErrorBanner.vue'
import { useMoleculeStore } from 'src/stores/molecule'

const molStore = useMoleculeStore()

</script>
