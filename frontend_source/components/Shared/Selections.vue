<template>
 <div class='q-mt-lg options'>

   <div class="col-6 items-center">

     <div v-for="item in props.store.options.slice(0, -1)">
       <q-checkbox
         dense
         class="q-mr-lg"
         :key="item.label"
         :val="item.value"
         :label="item.label"
         v-model="props.store.selection"
         color="blue"
       >

         <q-tooltip
           v-if="!props.store.hideOptionsTooltips"
           class="text-caption"
           anchor="bottom middle"
           :offset="[0, 1]"
           transition-show="scale"
           transition-hide="scale"
         >
           {{ item.tooltip }}
         </q-tooltip>
       </q-checkbox>
     </div>

   </div>

   <div class="col-6 items-center">
     <div v-for="item in props.store.options.slice(-1)">
       <q-checkbox
         dense
         class="q-mr-lg"
         :key="item.label"
         :val="item.value"
         :label="item.label"
         v-model="props.store.selection"
         color="blue"
       >

         <q-tooltip
           v-if="!props.store.hideOptionsTooltips"
           class="text-caption"
           anchor="bottom middle"
           :offset="[0, 1]"
           transition-show="scale"
           transition-hide="scale"
         >
           {{ item.tooltip }}
         </q-tooltip>
       </q-checkbox>
     </div>
     <div class="row items-center q-gutter-sm">
       <q-input
         dense
         stack-label
         v-model="props.store.angle"
         placeholder="0.0"
         type="number"
         label="rotate"
         style="max-width: 70px;"
       />
     </div>
     <div class="row items-center q-gutter-sm">
       <q-input
         label="indentation"
         placeholder="4"
         type="number"
         v-model="props.store.indentation"
         style="max-width: 70px;"
       />
     </div>
     <q-select
       stack-label
       dense
       options-dense
       hide-bottom-space
       hide-dropdown-icon
       transition-show="scale"
       transition-hide="scale"
       label="hydrogens"
       v-model="props.store.h2"
       :options="props.store.h2_options"
       style="width: 70px;"
     />

   </div>


  </div>
</template>

<script setup>

const props = defineProps([ "store" ])

</script>

<style scoped>
.options {
  display: flex;
  justify-content: space-evenly;
}
</style>
