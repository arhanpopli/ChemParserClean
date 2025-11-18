<template>
  <q-list
    bordered
    separator
    padding
    class="qlist"
  >
    <q-item-label
      header
      class='text-center text-h6 text-bold'
    > {{ props.header }}
    </q-item-label>

    <q-item
      tag="label"
      v-for="item in props.settings"
      :key="item.id_"
    >
      <q-item-section>
        <q-item-label
          class="text-subtitle1"
        >
          {{ item.description }}
        </q-item-label>
      </q-item-section>
      <q-item-section side >
        <q-toggle
          checked-icon="check"
          unchecked-icon="clear"
          color="blue"
          v-model="item.model"
          @update:model-value="OnChange(item)"
        />
      </q-item-section>
    </q-item>

  </q-list>
</template>

<script setup>

import { computed } from 'vue'

const props = defineProps(
  [
    "settings",
    "store",
    "header",
    "settingsStore"
  ]
)

const OnChange = (item) => {
  props.settingsStore.setValue(item)
  if (item['tag'] == "hideTooltips") {
    props.store.hideTooltips(item.model)
  }
}

</script>

<style scoped>
.qlist {
  width:70%;
  margin: 0 auto;
  margin-top: 30px;
}
</style>
