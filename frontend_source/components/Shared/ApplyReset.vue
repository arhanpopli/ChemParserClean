<template>
  <div class='row q-mt-md justify-center q-gutter-sm'>
    <Button
      color="positive"
      label="Apply"
      @click="props.apply"
      :loading="props.isApplying"
    />

    <Button
      color="negative"
      label="Reset"
      @click="props.reset"
      :loading="props.isResetting"
    />
  </div>
</template>

<script setup>

import Button from './Button.vue'

const props = defineProps({
  apply: {
    type: Function,
    required: true
  },
  reset: {
    type: Function,
    required: true
  },
  isApplying: {
    type: Boolean,
    required: true
  },
  isResetting: {
    type: Boolean,
    required: true
  }

})
</script>
