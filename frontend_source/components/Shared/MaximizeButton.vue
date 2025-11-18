<template>
  <div class="maxButton">
    <q-btn
      round
      flat
      class="text-grey-6"
      icon="fullscreen"
      size="lg"
      @click="props.store.maximizeComposer"
      style="z-index: 2;position: absolute"
    >
      <q-tooltip
        anchor="top middle"
        :offset="[1, 30]"
        class="bg-grey-9 text-subtitle2"
      >
        Maximize composer
      </q-tooltip>
    </q-btn>
  </div>
</template>

<script setup>


const props = defineProps(['store'])

</script>

<style>
.maxButton {
  margin-bottom: 10px;
  margin-left: 93%;
  display: flex;
}
</style>
