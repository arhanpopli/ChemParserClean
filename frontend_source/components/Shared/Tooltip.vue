<template>
  <div class="float-right q-mt-lg">
    <q-tooltip
      v-model="props.store.showCopyTooltip"
      :offset="[100, 20]"

    >
      Copied!</q-tooltip>
  </div>
</template>

<script setup>


const props = defineProps(['store'])

</script>

<style scoped>

</style>
