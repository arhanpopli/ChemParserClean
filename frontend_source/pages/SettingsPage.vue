<template>
  <q-page padding>
    <SettingsList
      :settings="molSettings"
      :store="molStore"
      :header="molHeader"
      :settingsStore="settingsStore"
    />

    <SettingsList
      :settings="reactionSettings"
      :store="reactionStore"
      :header="reactionHeader"
      :settingsStore="settingsStore"
    />
  </q-page>
</template>

<script setup>

import SettingsList from 'components/Settings/SettingsList.vue'
import { useSettings } from 'src/stores/settings'
import { useMoleculeStore } from 'src/stores/molecule'
import { useReactionStore } from 'src/stores/reaction'

const settingsStore = useSettings()
const molStore = useMoleculeStore()
const reactionStore = useReactionStore()

// Molecule's Settings section
const molSettings = settingsStore.settings.filter(item => item['type'] == "molecule")
const molHeader = "Molecule's Settings"

// Reaction's Settings section
const reactionSettings = settingsStore.settings.filter(item => item['type'] == "reaction")
const reactionHeader = "Reaction's Settings"

</script>
