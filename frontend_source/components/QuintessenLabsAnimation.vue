<template>
  <div class="quintessen-animation-wrapper">
    <div class="animation-container" :key="animationKey">
      <!-- Beaker SVGs -->
      <svg class="beaker to-q" viewBox="0 0 100 150" xmlns="http://www.w3.org/2000/svg">
        <defs>
          <linearGradient id="beakerGradient1" x1="0%" y1="0%" x2="0%" y2="100%">
            <stop offset="0%" style="stop-color:#e0e7ff;stop-opacity:0.8" />
            <stop offset="100%" style="stop-color:#c7d2fe;stop-opacity:0.9" />
          </linearGradient>
        </defs>
        <!-- Beaker body -->
        <path d="M 30 20 L 30 60 L 20 120 L 80 120 L 70 60 L 70 20 Z"
              fill="url(#beakerGradient1)"
              stroke="#a855f7"
              stroke-width="2"/>
        <!-- Beaker neck -->
        <rect x="35" y="10" width="30" height="10"
              fill="url(#beakerGradient1)"
              stroke="#a855f7"
              stroke-width="2"/>
        <!-- Liquid -->
        <path class="beaker-liquid" d="M 30 80 L 25 120 L 75 120 L 70 80 Z" />
      </svg>

      <svg class="beaker to-labs" viewBox="0 0 100 150" xmlns="http://www.w3.org/2000/svg">
        <defs>
          <linearGradient id="beakerGradient2" x1="0%" y1="0%" x2="0%" y2="100%">
            <stop offset="0%" style="stop-color:#e0e7ff;stop-opacity:0.8" />
            <stop offset="100%" style="stop-color:#c7d2fe;stop-opacity:0.9" />
          </linearGradient>
        </defs>
        <path d="M 30 20 L 30 60 L 20 120 L 80 120 L 70 60 L 70 20 Z"
              fill="url(#beakerGradient2)"
              stroke="#64748b"
              stroke-width="2"/>
        <rect x="35" y="10" width="30" height="10"
              fill="url(#beakerGradient2)"
              stroke="#64748b"
              stroke-width="2"/>
        <path class="beaker-liquid" d="M 30 80 L 25 120 L 75 120 L 70 80 Z"
              style="fill: #64748b; opacity: 0.5"/>
      </svg>

      <svg class="beaker to-labs" viewBox="0 0 100 150" xmlns="http://www.w3.org/2000/svg">
        <defs>
          <linearGradient id="beakerGradient3" x1="0%" y1="0%" x2="0%" y2="100%">
            <stop offset="0%" style="stop-color:#e0e7ff;stop-opacity:0.8" />
            <stop offset="100%" style="stop-color:#c7d2fe;stop-opacity:0.9" />
          </linearGradient>
        </defs>
        <path d="M 30 20 L 30 60 L 20 120 L 80 120 L 70 60 L 70 20 Z"
              fill="url(#beakerGradient3)"
              stroke="#64748b"
              stroke-width="2"/>
        <rect x="35" y="10" width="30" height="10"
              fill="url(#beakerGradient3)"
              stroke="#64748b"
              stroke-width="2"/>
        <path class="beaker-liquid" d="M 30 80 L 25 120 L 75 120 L 70 80 Z"
              style="fill: #64748b; opacity: 0.5"/>
      </svg>

      <svg class="beaker to-labs" viewBox="0 0 100 150" xmlns="http://www.w3.org/2000/svg">
        <defs>
          <linearGradient id="beakerGradient4" x1="0%" y1="0%" x2="0%" y2="100%">
            <stop offset="0%" style="stop-color:#e0e7ff;stop-opacity:0.8" />
            <stop offset="100%" style="stop-color:#c7d2fe;stop-opacity:0.9" />
          </linearGradient>
        </defs>
        <path d="M 30 20 L 30 60 L 20 120 L 80 120 L 70 60 L 70 20 Z"
              fill="url(#beakerGradient4)"
              stroke="#64748b"
              stroke-width="2"/>
        <rect x="35" y="10" width="30" height="10"
              fill="url(#beakerGradient4)"
              stroke="#64748b"
              stroke-width="2"/>
        <path class="beaker-liquid" d="M 30 80 L 25 120 L 75 120 L 70 80 Z"
              style="fill: #64748b; opacity: 0.5"/>
      </svg>

      <!-- Text elements -->
      <div class="letter-q">Q</div>
      <div class="quintessen-text">uintessen</div>
      <div class="labs-text">LABS</div>
    </div>

    <q-btn
      v-if="showRestartButton"
      class="restart-btn"
      label="Restart Animation"
      @click="restartAnimation"
      color="primary"
      rounded
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
    default: 8000 // 8 seconds
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
  min-height: 300px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  position: relative;
}

.animation-container {
  position: relative;
  width: 600px;
  height: 300px;
}

/* Beaker SVG Styles */
.beaker {
  position: absolute;
  width: 80px;
  height: 120px;
  opacity: 0;
  animation: fadeIn 0.6s ease-in-out forwards;
}

.beaker:nth-child(1) {
  left: 50px;
  top: 50%;
  transform: translateY(-50%);
  animation-delay: 0.2s;
}

.beaker:nth-child(2) {
  left: 150px;
  top: 50%;
  transform: translateY(-50%);
  animation-delay: 0.4s;
}

.beaker:nth-child(3) {
  left: 250px;
  top: 50%;
  transform: translateY(-50%);
  animation-delay: 0.6s;
}

.beaker:nth-child(4) {
  left: 350px;
  top: 50%;
  transform: translateY(-50%);
  animation-delay: 0.8s;
}

/* Beaker to Q transformation */
.beaker.to-q {
  animation: fadeIn 0.6s ease-in-out forwards, morphToQ 1s ease-in-out 1.5s forwards;
}

/* Beakers to "Labs" transformation */
.beaker.to-labs {
  animation: fadeIn 0.6s ease-in-out forwards, morphToLabs 1s ease-in-out 1.5s forwards;
}

/* The Q letter */
.letter-q {
  position: absolute;
  left: 50px;
  top: 50%;
  transform: translateY(-50%);
  font-size: 80px;
  font-weight: bold;
  color: #a855f7;
  opacity: 0;
  animation: revealQ 0.8s ease-in-out 2.5s forwards;
}

/* The "uintessen" text */
.quintessen-text {
  position: absolute;
  left: 140px;
  top: 50%;
  transform: translateY(-50%);
  font-size: 48px;
  font-weight: 600;
  color: #a855f7;
  opacity: 0;
  animation: slideOutFromQ 1s ease-out 3s forwards;
}

/* The "Labs" text */
.labs-text {
  position: absolute;
  left: 50px;
  top: calc(50% + 60px);
  font-size: 48px;
  font-weight: 300;
  color: #64748b;
  opacity: 0;
  letter-spacing: 8px;
  animation: revealLabs 0.8s ease-in-out 2.5s forwards;
}

/* Liquid in beakers */
.beaker-liquid {
  fill: #a855f7;
  opacity: 0.7;
  animation: bubble 2s ease-in-out infinite;
}

/* Animations */
@keyframes fadeIn {
  from {
    opacity: 0;
    transform: translateY(-50%) scale(0.5);
  }
  to {
    opacity: 1;
    transform: translateY(-50%) scale(1);
  }
}

@keyframes morphToQ {
  0% {
    opacity: 1;
    transform: translateY(-50%) scale(1);
  }
  50% {
    opacity: 0.3;
    transform: translateY(-50%) scale(1.2) rotate(180deg);
  }
  100% {
    opacity: 0;
    transform: translateY(-50%) scale(0.1) rotate(360deg);
  }
}

@keyframes morphToLabs {
  0% {
    opacity: 1;
    transform: translateY(-50%) scale(1);
  }
  50% {
    opacity: 0.3;
    transform: translateY(-50%) scale(1.2);
  }
  100% {
    opacity: 0;
    transform: translateY(-50%) translateX(50px) scale(0.1);
  }
}

@keyframes revealQ {
  from {
    opacity: 0;
    transform: translateY(-50%) scale(0.5) rotate(-90deg);
  }
  to {
    opacity: 1;
    transform: translateY(-50%) scale(1) rotate(0deg);
  }
}

@keyframes slideOutFromQ {
  from {
    opacity: 0;
    transform: translateY(-50%) translateX(-80px);
    letter-spacing: -10px;
  }
  to {
    opacity: 1;
    transform: translateY(-50%) translateX(0);
    letter-spacing: 2px;
  }
}

@keyframes revealLabs {
  from {
    opacity: 0;
    transform: translateY(20px);
    letter-spacing: 20px;
  }
  to {
    opacity: 1;
    transform: translateY(0);
    letter-spacing: 8px;
  }
}

@keyframes bubble {
  0%, 100% {
    transform: translateY(0);
  }
  50% {
    transform: translateY(-3px);
  }
}

/* Restart button */
.restart-btn {
  margin-top: 40px;
}

/* Glow effect */
.letter-q, .quintessen-text {
  text-shadow: 0 0 20px rgba(168, 85, 247, 0.5);
}

.labs-text {
  text-shadow: 0 0 15px rgba(100, 116, 139, 0.3);
}

/* Responsive design */
@media (max-width: 768px) {
  .animation-container {
    width: 100%;
    max-width: 400px;
    height: 200px;
    transform: scale(0.7);
  }

  .letter-q {
    font-size: 60px;
  }

  .quintessen-text {
    font-size: 36px;
  }

  .labs-text {
    font-size: 36px;
  }
}
</style>
