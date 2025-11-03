<template>
  <div>
    <q-btn-toggle
      v-model="model"
      spread
      size="12px"
      toggle-color="orange"
      :options="navLinks"
      :to="navLinks.path"
    />
  </div>
</template>

<script setup>

import { ref } from 'vue'

const navLinks = [
  {
    label: "Molecule",
    value: 'molecule',
    path: '/'
  },
  {
    label: "Reaction",
    value: 'reaction',
    path: "/mol_2_chemfig/reactions"
  },
  {
    label: "Tutorial",
    value: "tutorial",
    path: "/mol_2_chemfig/tutorial"
  },
  {
    label: 'About',
    value: 'about',
    path: '/mol_2_chemfig/about'
  },
  {
    label: 'Settings',
    value: 'settings',
    path: '/mol_2_chemfig/settings'
  }
]

const model = ref('molecule')

</script>
