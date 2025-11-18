<template>
  <div class="row justify-center text-center">
    <q-banner
      v-if="props.store.error"
      rounded
      dense
      inline-actions
      class="text-red-9 bg-grey-4"
      style="width:40%;"
    >
      {{ props.store.error }}
      <template v-slot:action>
        <q-btn
          dense
          flat
          icon="close"
          color="red"
          @click="props.store.closeBanner"
        />
      </template>
    </q-banner>
  </div>
</template>

<script setup>
const props = defineProps(['store'])
</script>
