<template>
  <div class="quintessen-animation-wrapper">
    <div class="animation-container" :key="animationKey">
      <svg id="morphingSVG" viewBox="0 0 900 500" xmlns="http://www.w3.org/2000/svg">
        <!-- Round-bottom flask (morphs to Q) -->
        <g class="glass-1">
          <path d="M 150 80 L 150 120 M 140 120 L 160 120 M 150 120 Q 120 150 120 180 Q 120 220 150 240 Q 180 220 180 180 Q 180 150 150 120">
            <animate
              attributeName="d"
              dur="1.2s"
              begin="2s"
              fill="freeze"
              values="
                M 150 80 L 150 120 M 140 120 L 160 120 M 150 120 Q 120 150 120 180 Q 120 220 150 240 Q 180 220 180 180 Q 180 150 150 120;
                M 150 100 Q 130 100 120 120 Q 110 140 110 165 Q 110 190 120 210 Q 130 230 150 230 Q 170 230 180 210 Q 190 190 190 165 Q 190 140 180 120 Q 170 100 150 100;
                M 150 120 Q 120 120 120 160 Q 120 200 150 200 Q 180 200 180 160 Q 180 120 150 120 M 170 185 L 185 205
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2s" fill="freeze" values="1; 1; 1" />
          </path>
        </g>

        <!-- Erlenmeyer flask (morphs to part of LABS - L) -->
        <g class="glass-2">
          <path d="M 320 280 L 320 240 L 310 240 L 310 200 L 330 200 L 330 240 L 320 240 L 340 300 L 300 300 Z">
            <animate
              attributeName="d"
              dur="1.2s"
              begin="2.2s"
              fill="freeze"
              values="
                M 320 280 L 320 240 L 310 240 L 310 200 L 330 200 L 330 240 L 320 240 L 340 300 L 300 300 Z;
                M 320 270 L 320 250 L 315 240 L 325 240 L 320 250 L 335 290 L 305 290 Z;
                M 150 280 L 150 350 L 200 350
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2.2s" fill="freeze" values="1; 1; 1" />
          </path>
        </g>

        <!-- Flat-bottom beaker (morphs to part of LABS - A) -->
        <g class="glass-3">
          <path d="M 400 220 L 400 300 L 460 300 L 460 220 M 405 220 L 455 220">
            <animate
              attributeName="d"
              dur="1.2s"
              begin="2.4s"
              fill="freeze"
              values="
                M 400 220 L 400 300 L 460 300 L 460 220 M 405 220 L 455 220;
                M 405 230 L 405 290 L 455 290 L 455 230;
                M 220 350 L 245 280 L 270 350 M 230 320 L 260 320
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2.4s" fill="freeze" values="1; 1; 1" />
          </path>
        </g>

        <!-- Volumetric flask (morphs to part of LABS - B) -->
        <g class="glass-4">
          <path d="M 540 200 L 540 220 M 535 220 L 545 220 M 540 220 Q 520 230 520 250 Q 520 270 540 280 Q 560 270 560 250 Q 560 230 540 220 M 530 280 L 530 300 L 550 300 L 550 280">
            <animate
              attributeName="d"
              dur="1.2s"
              begin="2.6s"
              fill="freeze"
              values="
                M 540 200 L 540 220 M 535 220 L 545 220 M 540 220 Q 520 230 520 250 Q 520 270 540 280 Q 560 270 560 250 Q 560 230 540 220 M 530 280 L 530 300 L 550 300 L 550 280;
                M 540 210 L 540 225 Q 525 235 525 250 Q 525 265 540 275 Q 555 265 555 250 Q 555 235 540 225 M 535 275 L 535 295 L 545 295 L 545 275;
                M 290 280 L 290 350 M 290 280 Q 315 280 315 295 Q 315 310 290 310 M 290 310 Q 315 310 315 325 Q 315 350 290 350
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2.6s" fill="freeze" values="1; 1; 1" />
          </path>
        </g>

        <!-- Additional flask for S in LABS -->
        <g class="glass-2" style="animation-delay: 1s;">
          <path d="M 620 240 L 620 280 Q 600 290 600 300 L 640 300 Q 640 290 620 280" opacity="0">
            <animate
              attributeName="d"
              dur="1s"
              begin="2.8s"
              fill="freeze"
              values="
                M 620 240 L 620 280 Q 600 290 600 300 L 640 300 Q 640 290 620 280;
                M 360 280 Q 340 280 340 295 Q 340 305 350 305 Q 360 305 360 315 Q 360 335 350 350 Q 340 350 340 340
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2.8s" fill="freeze" values="0; 1" />
          </path>
        </g>

        <!-- Sans-serif Q (final form) -->
        <text x="110" y="210" font-family="Arial, Helvetica, sans-serif" font-size="120" font-weight="bold" class="text-path">
          Q
          <animate attributeName="opacity" dur="0.5s" begin="3.2s" fill="freeze" values="0; 1" />
        </text>

        <!-- "uintessen" slides out from Q -->
        <text x="230" y="210" font-family="Arial, Helvetica, sans-serif" font-size="90" font-weight="300" class="text-path" opacity="0">
          uintessen
          <animate attributeName="opacity" dur="0.8s" begin="3.5s" fill="freeze" values="0; 1" />
          <animateTransform
            attributeName="transform"
            type="translate"
            dur="0.8s"
            begin="3.5s"
            fill="freeze"
            values="-60,0; 0,0"
          />
        </text>

        <!-- LABS (final form) -->
        <text x="110" y="380" font-family="Arial, Helvetica, sans-serif" font-size="90" font-weight="300" letter-spacing="15" class="text-path">
          LABS
          <animate attributeName="opacity" dur="0.5s" begin="3.2s" fill="freeze" values="0; 1" />
        </text>
      </svg>
    </div>

    <q-btn
      v-if="showRestartButton"
      class="restart-btn"
      label="Restart Animation"
      @click="restartAnimation"
      flat
      color="dark"
    />
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted } from 'vue'

// Props
const props = defineProps({
  autoRestart: {
    type: Boolean,
    default: true
  },
  restartInterval: {
    type: Number,
    default: 8000
  },
  showRestartButton: {
    type: Boolean,
    default: true
  }
})

// State
const animationKey = ref(0)
let intervalId = null

// Methods
const restartAnimation = () => {
  animationKey.value++
}

// Lifecycle hooks
onMounted(() => {
  if (props.autoRestart) {
    intervalId = setInterval(() => {
      restartAnimation()
    }, props.restartInterval)
  }
})

onUnmounted(() => {
  if (intervalId) {
    clearInterval(intervalId)
  }
})

// Expose methods for parent component
defineExpose({
  restartAnimation
})
</script>

<style scoped>
.quintessen-animation-wrapper {
  width: 100%;
  min-height: 500px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  position: relative;
  background: #ffffff;
}

.animation-container {
  position: relative;
  width: 900px;
  height: 500px;
  max-width: 100%;
}

#morphingSVG {
  width: 100%;
  height: 100%;
}

/* Lab glassware - black stroke only */
.glass-1,
.glass-2,
.glass-3,
.glass-4 {
  fill: none;
  stroke: #000000;
  stroke-width: 2.5;
  stroke-linecap: round;
  stroke-linejoin: round;
  opacity: 0;
  animation: fadeInGlass 0.8s ease-out forwards;
}

.glass-1 { animation-delay: 0.2s; }
.glass-2 { animation-delay: 0.4s; }
.glass-3 { animation-delay: 0.6s; }
.glass-4 { animation-delay: 0.8s; }

/* Sans-serif text */
.text-path {
  fill: #000000;
  opacity: 0;
}

@keyframes fadeInGlass {
  to {
    opacity: 1;
  }
}

/* Restart button */
.restart-btn {
  margin-top: 40px;
  border: 1px solid #000000;
}

/* Responsive design */
@media (max-width: 768px) {
  .animation-container {
    width: 100%;
    max-width: 600px;
    height: 350px;
  }

  #morphingSVG {
    transform: scale(0.7);
  }
}

@media (max-width: 480px) {
  .animation-container {
    width: 100%;
    height: 300px;
  }

  #morphingSVG {
    transform: scale(0.5);
  }
}
</style>
