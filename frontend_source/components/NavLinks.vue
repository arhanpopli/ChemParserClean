<template>
  <div>
    <q-btn v-for="link in navLinks"
      :key="link.label"
      flat
      color="white"
      :label="link.label"
      :to="link.path"
    />
  </div>
</template>

<script setup>

const navLinks = [
  {
    label: "Molecule",
    path: '/'
  },
  {
    label: "Reaction",
    path: "/reaction"
  },
  {
    label: "Tutorial",
    path: "/tutorial"
  },
  // {
  //   label: 'About',
  //   path: '/about'
  // },
  {
    label: 'Settings',
    path: '/settings'
  }
]

</script>
