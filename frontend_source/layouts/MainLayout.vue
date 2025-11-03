<template>
  <q-layout view="lHh Lpr lFf">
    <q-header elevated>
      <q-toolbar
          style="min-height:35px"
      >
        <q-toolbar-title
          class="text-bold q-ml-xl"
        >
          Mol2chemfig Web
        </q-toolbar-title>

        <div><NavLinks /></div>
      </q-toolbar>
    </q-header>

    <q-page-container>
      <router-view />
    </q-page-container>

    <q-footer elevated>
      <q-toolbar
        style="min-height:30px"
      >
        <q-space />
        Mol2Chemfig Web &copy Vitali Coltuclu <span class="q-ml-md">2014 - {{ currentYear }}</span>
        <q-space />
      </q-toolbar>
    </q-footer>
  </q-layout>
</template>


<script setup>
import { computed } from 'vue'
import NavLinks from 'components/NavLinks.vue'

const currentYear = computed(() => {
  const d = new Date();
  return d.getFullYear();
})

</script>
