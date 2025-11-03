<template>
  <router-view />
</template>

<script setup>
import { onBeforeMount } from 'vue'
import { useSettings } from 'src/stores/settings'

const settings = useSettings()


onBeforeMount(() => {
      settings.setSettings()
    })
</script>
